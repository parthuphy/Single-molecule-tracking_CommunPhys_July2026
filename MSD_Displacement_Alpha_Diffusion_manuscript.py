# =============================================================================
# Purpose:
# This code calculates the mean squared displacement (MSD) of individual tracked
# particles from single-molecule tracking (SMT) trajectory data. It also
# computes the ensemble-averaged MSD, particle displacement statistics,
# anomalous diffusion exponent (α), diffusion coefficient (D), and summary
# statistics (mean, standard deviation, and standard error of D).

# Required input:
# A CSV file containing the trajectory data with the following header columns:
#
#     TRACK_ID, POSITION_X, POSITION_Y, FRAME
#
# where:
#     TRACK_ID   : Unique particle ID
#     POSITION_X : X-coordinate of the particle (µm)
#     POSITION_Y : Y-coordinate of the particle (µm)
#     FRAME      : Frame number
# =============================================================================

import pandas as pd
import numpy as np

frame_rate = 10.0  # Frames per second
min_particles_for_avg = 1  # Minimum number of particles for ensemble MSD at a lag time
tau_cutoff = 1  # The straight line will be fitted upto tau_cutoff seconds for calculating the exponent in log-log plot of MSD vs. lag time

# Provide the path of input csv file
input_file = r'input_csv_file_path'

# Provide the path of output excel file
output_file = r'output_excel_file_path'

# Load the input data
df = pd.read_csv(input_file, skiprows=1) # Skipping first row of the input csv file as it contains header only
df.columns = ['TRACK_ID', 'X', 'Y', 'FRAME']

# Calculation of MSD and Displacement of individual tracked particles
msd_data = {}
displacement_records = []

for track_id, group in df.groupby('TRACK_ID'):
    group = group.sort_values('FRAME')
    coords = group[['X', 'Y']].to_numpy()
    frames = group['FRAME'].values
    n = len(coords) # Number of frames for the current tracked particle

    lags = np.arange(1, n)
    lag_times_sec = lags / frame_rate
    msds = []

    for lag in lags:
        displacements = coords[lag:] - coords[:-lag] # Displacement calculation
        squared_displacement = np.sum(displacements**2, axis=1)
        msd = np.mean(squared_displacement) # MSD calculation
        msds.append(msd)

    # Summary table for MSD vs. lag time
    msd_data[track_id] = pd.DataFrame({
        f'Lag_{track_id} (s)': lag_times_sec,
        f'MSD_{track_id}': msds
    })

    # Δx, Δy, Δr between consecutive frames
    delta_coords = coords[1:] - coords[:-1]
    delta_r = np.linalg.norm(delta_coords, axis=1)
    delta_x = delta_coords[:, 0]
    delta_y = delta_coords[:, 1]

    for i in range(len(delta_r)):
        displacement_records.append({
            'TRACK_ID': track_id,
            'Time_i (s)': frames[i] / frame_rate,
            'Time_i+1 (s)': frames[i + 1] / frame_rate,
            'Δx': delta_x[i],
            'Δy': delta_y[i],
            'Δr': delta_r[i]
        })

# Combine MSDs
df_msd = pd.concat(msd_data.values(), axis=1)
df_disp = pd.DataFrame(displacement_records)

# Count number of tracked particles
num_particles_tracked = len(msd_data)
summary_df = pd.DataFrame({'Particles with MSD': [num_particles_tracked]})

# Calculation of Ensemble-averaged MSD and number of tracked particles for a given lag time
max_lag_len = max([df_msd.iloc[:, i].dropna().shape[0]
                   for i in range(0, df_msd.shape[1], 2)])

ensemble_lags = []
ensemble_msds = []
ensemble_counts = []

for lag_idx in range(max_lag_len):
    msd_vals = []
    lag_val = None

    for i in range(0, df_msd.shape[1], 2):
        lag_col = df_msd.iloc[:, i]
        msd_col = df_msd.iloc[:, i + 1]

        # Check whether the current particle has a valid MSD value at the selected lag-time index
        if lag_idx < len(lag_col) and not pd.isna(msd_col.iloc[lag_idx]):
            msd_vals.append(msd_col.iloc[lag_idx])
            if lag_val is None:
                lag_val = lag_col.iloc[lag_idx]

    num_contrib_particles = len(msd_vals) # Number of particles contributing to this lag time

    # Calculate ensemble-averaged MSD only if sufficient number of particles are present in the system
    if num_contrib_particles >= min_particles_for_avg:
        avg_msd = np.mean(msd_vals) # Calculate ensemble-averaged MSD
        ensemble_lags.append(lag_val)
        ensemble_msds.append(avg_msd)
        ensemble_counts.append(num_contrib_particles)
    else:
        break

# Summary table for Ensemble-averaged MSD, lag time, and number of tracked particles
ensemble_df = pd.DataFrame({
    'Lag Time (s)': ensemble_lags,
    'Ensemble MSD': ensemble_msds,
    'N (particles)': ensemble_counts
})

# Separate N vs lag time sheet
N_vs_tau_df = ensemble_df[['Lag Time (s)', 'N (particles)']]

# Calculation of exponent values and Diffusion coefficient for each trajectories
alpha_records = []
D_values = []

for i in range(0, df_msd.shape[1], 2):

    lag_col = df_msd.iloc[:, i].dropna().values
    msd_col = df_msd.iloc[:, i + 1].dropna().values

    # Select only valid data points within the specified lag-time range and with positive MSD values
    valid_indices = (lag_col > 0) & (lag_col <= tau_cutoff) & (msd_col > 0)

    # Extract lag time and MSD values used for fitting
    lag_fit = lag_col[valid_indices]
    msd_fit = msd_col[valid_indices]

    # Extract particle (TRACK_ID) from the MSD column name
    msd_col_name = df_msd.columns[i + 1]
    track_id = msd_col_name.replace("MSD_", "")

    # Initialize fitting parameters
    alpha = np.nan
    intercept = np.nan
    D = np.nan

    if len(lag_fit) > 1: # Perform fitting only if at least two data points are available

        try:
            
            # Convert lag time and MSD to logarithmic scale
            log_tau = np.log10(lag_fit)
            log_msd = np.log10(msd_fit)

            # Perform linear fitting to log(MSD) = α log(τ) + log(4D) for 2D system
            slope, intercept = np.polyfit(log_tau, log_msd, 1)

            if slope >= 0: # Accept only physically meaningful (non-negative) exponent values

                alpha = slope

                # Diffusion coefficient
                D = (10 ** intercept) / 4 # Because Intercept = log10(4D)

                D_values.append(D)

        except:
            pass

    # Summary table for exponent values, particle ID, intercept and diffusion of coefficient values
    alpha_records.append({
        'Particle': track_id,
        'Alpha': alpha,
        'Intercept log10(4D)': intercept,
        'Diffusion Coefficient D (µm²/s)': D
    })

alpha_df = pd.DataFrame(alpha_records)

# Remove rows where alpha is NaN
alpha_df = alpha_df.dropna(subset=['Alpha'])

N = len(D_values) # Total number of particles with valid diffusion coefficient

D_mean = np.mean(D_values) # Mean diffusion coefficient

D_std = np.std(D_values, ddof=1) # Standard deviation of diffusion coefficient

D_se = D_std / np.sqrt(N) # Standard error of diffusion coefficient

# Summary table for Diffusion coefficient
D_summary_df = pd.DataFrame({
    'Quantity': [
        'Number of particles (N)',
        'Mean D',
        'Standard Deviation',
        'Standard Error'
    ],
    'Value': [
        N,
        D_mean,
        D_std,
        D_se
    ]
})

# Save the results
with pd.ExcelWriter(output_file) as writer:
    df_msd.to_excel(writer, index=False, sheet_name='MSD')
    df_disp.to_excel(writer, index=False, sheet_name='Displacements')
    ensemble_df.to_excel(writer, index=False, sheet_name='MSD_Average')
    alpha_df.to_excel(writer, index=False, sheet_name='Alpha_D_Values')
    D_summary_df.to_excel(writer, index=False, sheet_name='D_Statistics')
    summary_df.to_excel(writer, index=False, sheet_name='Summary')
    N_vs_tau_df.to_excel(writer, index=False, sheet_name='N_vs_tau')