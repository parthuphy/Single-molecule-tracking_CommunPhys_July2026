# =============================================================================
# Purpose:
# This code calculates the instantaneous velocity of particles from
# single-molecule tracking (SMT) data, filters particle positions based on a
# user-defined velocity threshold, and constructs a Voronoi tessellation using
# the retained particle positions. The code computes the area, perimeter, and
# number of neighboring cells for all bounded Voronoi polygons, exports these
# quantities to an Excel file, and generates both a full-field Voronoi diagram
# and a zoomed Voronoi diagram for a user-selected trajectory.
#
# Required Input:
# A CSV file containing the following columns:
#
#     TRACK_ID    : Unique identifier for each tracked particle
#     POSITION_X  : X-coordinate of the particle (µm)
#     POSITION_Y  : Y-coordinate of the particle (µm)
#     FRAME       : Frame number
# =============================================================================


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.spatial import Voronoi
from shapely.geometry import LineString, Polygon, box

# Input csv file path
csv_file = r'input_csv_file_path'

# Output excel file path
output_excel = r'output_excel_file_path'

dt = 0.1 # Seconds per frame (e.g, 10 frames = 1 second)
frame_step = 5 # Number of frames over which velocity will be calculated
velocity_threshold = 0.19823 # Maximum velocity of the particle (in µm/s) for which Voronoi tessellation will be constructed

trajectory_color = "blue" # Color of the trajectory segments
voronoi_color = "orange" # Color of the Voronoi tessellation lines

# Load input data
df = pd.read_csv(csv_file)
df = df.sort_values(by=["TRACK_ID", "FRAME"])

# Compute velocity over a given number of frame difference
df["dx"] = df.groupby("TRACK_ID")["POSITION_X"].diff(periods=frame_step)
df["dy"] = df.groupby("TRACK_ID")["POSITION_Y"].diff(periods=frame_step)
df["dframe"] = df.groupby("TRACK_ID")["FRAME"].diff(periods=frame_step)

df["dr"] = np.sqrt(df["dx"]**2 + df["dy"]**2) # Calculation of the square of displacement
df["dt"] = df["dframe"] * dt # Convert frame number into time (s)
df["velocity"] = df["dr"] / df["dt"] # Calculate instantaneous velocity

# Assign velocity to the earlier frame of the particle
df["velocity"] = df.groupby("TRACK_ID")["velocity"].shift(-frame_step)

# Keep only trajectries with selected velocity range
df = df[np.isfinite(df["velocity"])]
df = df[df["velocity"] > 0]
df = df[df["velocity"] <= velocity_threshold]

# Plot filtered particle trajectories
fig, ax = plt.subplots(figsize=(8, 6)) # Figure for plotting trajectories and Voronoi tessellation

for _, track in df.groupby("TRACK_ID"):
    x = track["POSITION_X"].values
    y = track["POSITION_Y"].values

    # Skip trajectories having fewer than two points
    if len(x) < 2:
        continue

    # Convert trajectory coordinates into connected line segments
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, colors=trajectory_color, linewidths=3, zorder=3) # Draw the trajectory
    ax.add_collection(lc)

# Construction of Voronoi tessellation from filtered trajectory points
# Extract particle coordinates used as Voronoi seed points
points_xy = df[["POSITION_X", "POSITION_Y"]].to_numpy()

# Determine the spatial limits of the trajectory data
xmin, xmax = points_xy[:, 0].min(), points_xy[:, 0].max()
ymin, ymax = points_xy[:, 1].min(), points_xy[:, 1].max()

# Define padding around the data to create a clipping boundary
padding = 1.0
bbox = box(xmin - padding, ymin - padding,
           xmax + padding, ymax + padding)

cell_areas = []
cell_perimeters = []
neighbor_counts = []

# Construct Voronoi tessellation only if sufficient seed points are available
if len(points_xy) >= 4:

    vor = Voronoi(points_xy)

    # Count the number of neighboring cells sharing a Voronoi edge
    neighbor_dict = {i: 0 for i in range(len(points_xy))}
    for p1, p2 in vor.ridge_points:
        neighbor_dict[p1] += 1
        neighbor_dict[p2] += 1

    for point_idx, region_idx in enumerate(vor.point_region):

        region = vor.regions[region_idx]

        # Remove unbounded cells
        if -1 in region or len(region) == 0:
            continue

        # Construct the polygon corresponding to the Voronoi cell
        polygon_coords = [vor.vertices[i] for i in region]
        poly = Polygon(polygon_coords)

        # Ignore invalid polygons
        if not poly.is_valid:
            continue

        # Store area, perimeter and number of neighbors
        cell_areas.append(poly.area)
        cell_perimeters.append(poly.length)
        neighbor_counts.append(neighbor_dict[point_idx])

    # Plot Voronoi edges
    for ridge_vertices in vor.ridge_vertices:

        # Ignore infinite Voronoi edges
        if -1 in ridge_vertices:
            continue

        v0, v1 = vor.vertices[ridge_vertices]

        # Construct the Voronoi edge
        line = LineString([v0, v1])
        clipped = line.intersection(bbox) # Clip the edge to the plotting boundary

        # Ignore empty intersections
        if clipped.is_empty:
            continue

        # Plot the clipped Voronoi edge
        if clipped.geom_type == "LineString":
            x_clip, y_clip = clipped.xy
            ax.plot(x_clip, y_clip,
                    linestyle="--",
                    color=voronoi_color,
                    linewidth=1.5,
                    alpha=1,
                    zorder=1)

# Source data for Voronoi figure
voronoi_seed_points = (
    df[["TRACK_ID", "FRAME", "POSITION_X", "POSITION_Y", "velocity"]]
    .rename(columns={
        "POSITION_X": "X (µm)",
        "POSITION_Y": "Y (µm)",
        "velocity": "Velocity (µm/s)"
    })
    .sort_values(["TRACK_ID", "FRAME"])
    .reset_index(drop=True)
)

# Export Voronoi data to an excel file
with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:

    # Source data used to construct the figure
    voronoi_seed_points.to_excel(
        writer,
        sheet_name="Voronoi_Seed_Points",
        index=False
    )

    # Voronoi Cell area
    pd.DataFrame({
        "Cell_Area_um2": cell_areas
    }).to_excel(
        writer,
        sheet_name="Cell_Area",
        index=False
    )

    # Voronoi Cell perimeters
    pd.DataFrame({
        "Cell_Perimeter_um": cell_perimeters
    }).to_excel(
        writer,
        sheet_name="Cell_Perimeter",
        index=False
    )

    # Number of neighbors
    pd.DataFrame({
        "Number_of_Neighbors": neighbor_counts
    }).to_excel(
        writer,
        sheet_name="Num_Neighbors",
        index=False
    )

# Final plot formatting
ax.set_xlabel("X [µm]") # Label X-axis
ax.set_ylabel("Y [µm]") # Label Y-axis
ax.set_title("2D trajectories with Voronoi tessellation") # Title of the Voronoi plot
ax.set_aspect("equal") # Keep identical scaling on both axis

# Set plotting limits
ax.set_xlim(xmin - padding, xmax + padding)
ax.set_ylim(ymin - padding, ymax + padding)

plt.tight_layout()

# Save the figure with high resolutiion
plt.savefig(
    r'output_png_file_path',
    dpi=1200,
    bbox_inches="tight"
)
plt.show()

# Zoomed trajectory and Voronoi diagram

# Rank all trajectories according to the number of recorded positions
trajectory_lengths = (
    df.groupby("TRACK_ID")
      .size()
      .sort_values(ascending=False)
)

ranked_tracks = trajectory_lengths.reset_index()
ranked_tracks.columns = ["TRACK_ID", "Length"]

# Assign a rank to each trajectory (1 = longest trajectory)
ranked_tracks["Rank"] = np.arange(1, len(ranked_tracks) + 1)

# Provide the trajectory rank that need to be plotted
selected_rank = 2

# Check that the selected rank is valid
if selected_rank < 1 or selected_rank > len(ranked_tracks):
    raise ValueError("Invalid rank selected.")

# Obtain the TRACK_ID corresponding to the selected rank
selected_track_id = ranked_tracks.loc[
    ranked_tracks["Rank"] == selected_rank, "TRACK_ID"
].values[0]

# Extract selected trajectory
zoom_df = df[df["TRACK_ID"] == selected_track_id]

# Extract trajectory coordinates for Voronoi construction
points_zoom = zoom_df[["POSITION_X", "POSITION_Y"]].to_numpy()

# Plot the selected trajectory

# Create the figure of zoomed trajectory
fig_zoom, ax_zoom = plt.subplots(figsize=(6, 6))

# Extract X- and Y-coordinates of the selected trajectory
x = zoom_df["POSITION_X"].values
y = zoom_df["POSITION_Y"].values

# Plot the trajectory only if at least two positions are available
if len(x) >= 2:

    # Convert trajectory coordinates into connected line segments
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Draw the trajectories
    lc = LineCollection(segments, colors="blue", linewidths=2, alpha=0.7, zorder=3)
    ax_zoom.add_collection(lc)

# Plot individual trajectory points
ax_zoom.scatter(x, y, color="black", s=2, zorder=4)

# Construct Voronoi tesselation from the selected trajectory

if len(points_zoom) >= 4: # Construct the Voronoi tessellation only if sufficient points are available

    vor_zoom = Voronoi(points_zoom)

    # Determine the spatial limits of the selected trajectory
    xmin_z, xmax_z = points_zoom[:, 0].min(), points_zoom[:, 0].max()
    ymin_z, ymax_z = points_zoom[:, 1].min(), points_zoom[:, 1].max()

    # Calculate the spatial extent of the trajectory
    x_span = xmax_z - xmin_z
    y_span = ymax_z - ymin_z

    # Define a small padding around the trajectory
    padding_x = 0.05 * x_span if x_span > 0 else 0.01
    padding_y = 0.05 * y_span if y_span > 0 else 0.01

    # Construct a tight bounding box for clipping Voronoi edges
    bbox_z = box(xmin_z - padding_x,
                ymin_z - padding_y,
                xmax_z + padding_x,
                ymax_z + padding_y)

    # Identify bounded Voronoi regions only
    bounded_regions = set()

    for point_idx, region_idx in enumerate(vor_zoom.point_region):
        region = vor_zoom.regions[region_idx]

        # Retain only finite Voronoi cells
        if -1 not in region and len(region) > 0:
            bounded_regions.add(region_idx)

    for (p1, p2), ridge_vertices in zip(vor_zoom.ridge_points,
                                        vor_zoom.ridge_vertices):

        # Ignore infinite Voronoi edges
        if -1 in ridge_vertices:
            continue

        # Identify the Voronoi regions sharing the ridge
        r1 = vor_zoom.point_region[p1]
        r2 = vor_zoom.point_region[p2]

        # Plot only edges shared by two bounded Voronoi cells
        if r1 not in bounded_regions or r2 not in bounded_regions:
            continue

        # Construct the Voronoi edge
        v0, v1 = vor_zoom.vertices[ridge_vertices]
        line = LineString([v0, v1])

        # Clip the edge to the plotting boundary
        clipped = line.intersection(bbox_z)

        # Ignore empty intersections
        if clipped.is_empty:
            continue

        # Plot the clipped Voronoi edge
        if clipped.geom_type == "LineString":
            x_clip, y_clip = clipped.xy
            ax_zoom.plot(
                x_clip,
                y_clip,
                linestyle="--",
                color="orange",
                linewidth=1,
                alpha=1,
                zorder=1)

# Final formatting of zoomed trajectory figure
ax_zoom.set_xlabel("X [µm]") # Label X-axis
ax_zoom.set_ylabel("Y [µm]") # Label Y-axis
ax_zoom.set_title(f"Trajectory with (Rank {selected_rank}) and Voronoi") # Title of the figure
ax_zoom.set_aspect("equal")

# Set the plotting limits
ax_zoom.set_xlim(xmin_z - padding_x, xmax_z + padding_x)
ax_zoom.set_ylim(ymin_z - padding_y, ymax_z + padding_y)

plt.tight_layout()

# Output file path
zoom_output_path = r'output_png_zoomed_file_path'

# Save the high resolution figure
plt.savefig(zoom_output_path, dpi=1200, bbox_inches="tight")

# Display the figure
plt.show()