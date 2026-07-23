# =============================================================================
# Purpose:
# This code calculates the two-point overlap correlation function, Q₂(a,t),
# from single-molecule tracking (SMT) trajectory data over a user-defined range
# of probe lengths (a). The calculated Q₂ values as a function of delay time
# are exported to an Excel workbook, with each selected probe length saved in
# a separate worksheet.
#
# Required input:
# A CSV file containing the trajectory data with the following header columns:
#
#     TRACK_ID, POSITION_X, POSITION_Y, FRAME
#
# where:
#     TRACK_ID   : Unique particle/trajectory identifier
#     POSITION_X : X-coordinate of the particle (µm)
#     POSITION_Y : Y-coordinate of the particle (µm)
#     FRAME      : Frame number
# =============================================================================

import pandas as pd
import numpy as np
import os

# Input csv file path
input_csv = r'input_csv_file_path'

# Output excel file path
output_excel = r'output_excel_file_path'

# Define the minimum and maximum probe length (a) in micrometers
a_min = 0.55
a_max = 0.65

delta_a = 0.05 # Step size used for generating probe length values

frame_rate = 10 # e.g., 10 frames = 1 second

# Load the input data
df = pd.read_csv(input_csv)

# Verify that the required columns are present in the input file
required_cols = ["TRACK_ID", "POSITION_X", "POSITION_Y", "FRAME"]
if not all(col in df.columns for col in required_cols):
    raise ValueError(f"Input CSV must contain columns: {required_cols}")

df["TIME"] = df["FRAME"] / frame_rate # Convert frame number to time (s)
df = df.sort_values(by=["TRACK_ID", "FRAME"]).reset_index(drop=True) # Sort trajectories by TRACK_ID and FRAME number
groups = df.groupby("TRACK_ID")

# Function to calculate the two-point overlap correlation function Q₂ for a given delay time and probe length (a)
def compute_Q2_for_delay(delay_t, a):
    Q_values = []

    for track_id, track in groups:
        times = track["TIME"].values
        x = track["POSITION_X"].values
        y = track["POSITION_Y"].values

        n_points = len(times) # Number of recorded positions for the current particle

        # Skip trajectories having fewer than two points
        if n_points < 2:
            continue
        
        for i, t0 in enumerate(times):
            t1 = t0 + delay_t # Calculate the target time (t₀ + delay time)
            j = np.where(np.isclose(times, t1, atol=1e-9))[0] # Find the corresponding frame at the target time
            
            # Skip if the target time is not present
            if len(j) == 0:
                continue

            j = j[0]

            # Calculate displacement components
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            dr2 = dx**2 + dy**2 

            Q_i = np.exp(-dr2 / (2 * a**2)) # Calculate Q₂ for the current displacement
            Q_values.append(Q_i)

    # Return NaN if no valid Q₂ values are found
    if len(Q_values) == 0:
        return np.nan

    return np.mean(Q_values)

# Computation of Q₂ for multiple probe lengths (a)
unique_times = sorted(df["TIME"].unique()) # Extract all time values
max_delay = unique_times[-1] - unique_times[0] # Maximum available delay time
delay_times = np.arange(0.1, max_delay + 0.1, 0.1) # Generate delay times with 0.1 s intervals

a_values = np.round(np.arange(a_min, a_max + delta_a / 2, delta_a), 3) # Generate probe length values within the specified range

results_dict = {}

for a in a_values:
    res = []
    for t in delay_times:
        Q2 = compute_Q2_for_delay(t, a) # Calculate Q₂
        res.append({"Delay_time_s": t, "Q2": Q2})
    df_res = pd.DataFrame(res)
    results_dict[a] = df_res

# User input for saving the data
save_all = input("\nDo you want to save ALL computed probe length (a) values? (yes/no): ").strip().lower()

if save_all in ["yes", "y"]:
    selected_a_values = list(results_dict.keys())
    print("\nSaving ALL computed probe lengths.")
else:
    user_input = input(
        "\nEnter the list of probe lengths (a) to save, separated by commas:\n"
        "Example: 0.20,0.25,0.30\n→ "
    )

    # Convert user input into numerical values
    try:
        selected_a_values = [round(float(x.strip()), 3) for x in user_input.split(",") if x.strip()]
    
    except ValueError:
        print("Invalid input. Please enter numeric values separated by commas.")
        exit()

    selected_a_values = [a for a in selected_a_values if a in results_dict]

    if len(selected_a_values) == 0:
        print("No valid probe lengths found in the computed range.")
        exit()

# Save selected data in excel file
os.makedirs(os.path.dirname(output_excel), exist_ok=True) # Create the output directory if it does not already exist

with pd.ExcelWriter(output_excel) as writer:
    for a in selected_a_values:
        results_dict[a].to_excel(writer, sheet_name=f"a={a:.3f}", index=False)