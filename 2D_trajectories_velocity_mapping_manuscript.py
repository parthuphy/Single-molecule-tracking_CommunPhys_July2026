# =============================================================================
# Purpose:
# This code calculates the instantaneous velocity of each tracked particle from
# its trajectory, computes the mean velocity, standard deviation,
# and standard error, and generates a 2D trajectory plot with velocity-based
# color mapping.
#
# Required input:
# A CSV file containing the trajectory data with the following header columns:
#     TRACK_ID, POSITION_X, POSITION_Y, FRAME
#
# where:
#     TRACK_ID   : Unique ID of each tracked particle
#     POSITION_X : X-coordinate (µm)
#     POSITION_Y : Y-coordinate (µm)
#     FRAME      : Frame number
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

# Provide the path of input csv file
csv_file = r'input_csv_file'

dt = 0.1   # Seconds per frame (e.g, 10 frames = 1 second)

frame_step = 5 # Number of frames over which velocity will be calculated

# Load input data
df = pd.read_csv(csv_file)

# Ensure correct sorting
df = df.sort_values(by=["TRACK_ID", "FRAME"])

# Differences within each track
df["dx"] = df.groupby("TRACK_ID")["POSITION_X"].diff(periods=frame_step)
df["dy"] = df.groupby("TRACK_ID")["POSITION_Y"].diff(periods=frame_step)
df["dframe"] = df.groupby("TRACK_ID")["FRAME"].diff(periods=frame_step)

# Displacement calculation
df["dr"] = np.sqrt(df["dx"]**2 + df["dy"]**2)

# Convert frame difference to time difference
df["dt"] = df["dframe"] * dt

# Velocity calculation
df["velocity"] = df["dr"] / df["dt"]

# Assign velocity to the earlier frame
df["velocity"] = (
    df.groupby("TRACK_ID")["velocity"]
      .shift(-frame_step)
)

# Remove NaN, -inf, +inf values of velocities
df = df[np.isfinite(df["velocity"])]

# Remove zero velocity values
df = df[df["velocity"] > 0]

# Remove invalid rows (first point of each track, gaps, etc.)
df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["velocity"])

# Mean velocity calculation
track_mean_velocity = df.groupby("TRACK_ID")["velocity"].mean()

# Number of tracked particles
N = len(track_mean_velocity)

# Mean, standard deviation, and standard error
mean_of_tracks = track_mean_velocity.mean()
std_of_tracks = track_mean_velocity.std(ddof=1)
se_of_tracks = std_of_tracks / np.sqrt(N)

# Print results
print(f"Number of tracked particles (N): {N}")
print(f"Mean velocity (track-averaged): {mean_of_tracks:.5f} µm/s")
print(f"Standard deviation: {std_of_tracks:.5f} µm/s")
print(f"Standard error: {se_of_tracks:.5f} µm/s")

# Plot trajectories with velocity color mapping
fig, ax = plt.subplots(figsize=(8, 6))

# Provide logarithmic range of velocity
norm = LogNorm(vmin=5e-3, vmax=5e0)

# Velocity color map
cmap = plt.cm.plasma

for track_id, track in df.groupby("TRACK_ID"):
    x = track["POSITION_X"].values
    y = track["POSITION_Y"].values
    v = track["velocity"].values

    # Skip trajectories having fewer than two points
    if len(x) < 2:
        continue

    # Construct line segments
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create a collection of coloured line segments
    lc = LineCollection(
        segments,
        cmap=cmap,
        norm=norm
    )
    # Assign velocity values to each line segment
    lc.set_array(v[:-1])

    # Set thickness of trajectory lines
    lc.set_linewidth(2)

    # Add trajectories to the plot
    ax.add_collection(lc)

# Add colorbar
cbar = plt.colorbar(lc, ax=ax)

# Label the colorbar
cbar.set_label("Velocity [µm/s]")

# Figure formatting
ax.set_xlabel("X [µm]") # X-axis label
ax.set_ylabel("Y [µm]") # Y-axis label
ax.set_title("2D-trajectories") # Title of the plot
ax.autoscale() 
ax.set_aspect("equal") # Keep identical scaling on X-Y axis

plt.tight_layout()

# Save figure with high resolution (provide the path of output png file)
plt.savefig(r'output_png_file', dpi=1200, bbox_inches="tight")

plt.show() # Display figure