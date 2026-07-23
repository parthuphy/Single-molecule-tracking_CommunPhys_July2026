# Single-molecule-tracking_CommunPhys_July2026
# This repository contains the Python scripts for single-molecule tracking analysis of thin polymer films.

> **Activated Cage-Escape Dynamics Bridge Molecular Fluctuations and Macroscopic Swelling in Glassy Polymer Films**

These scripts were developed to analyze molecular trajectories obtained from TrackMate and were used to calculate particle velocities, mean squared displacement (MSD), two-point overlap correlation functions, and Voronoi tessellation statistics.

---

## Repository Contents

| Script | Description |
|---------|-------------|
| `2D_trajectories_velocity_mapping_manuscript.py` | Calculates the instantaneous velocity of each tracked particle, computes the mean velocity, standard deviation, and standard error, and generates a 2D trajectory plot colored by particle velocity. |
| `MSD_Displacement_Alpha_Diffusion_manuscript.py` | Calculates individual and ensemble-averaged mean squared displacement (MSD), displacement statistics, anomalous diffusion exponent (α), diffusion coefficient (D), and exports all calculated quantities to an Excel workbook. |
| `Q2_calculation_manuscript.py` | Calculates the two-point overlap correlation function, Q₂(a,t), over a user-defined range of probe lengths and exports the results to an Excel workbook. |
| `2D_trajectory_Voronoi_manuscript.py` | Calculates instantaneous particle velocities, filters particle positions based on a user-defined velocity threshold, constructs Voronoi tessellations, computes bounded-cell statistics (area, perimeter, and number of neighbors), and generates full-field and zoomed Voronoi diagrams. |

---

# Input Data

All scripts require trajectory data exported from **TrackMate (FIJI/ImageJ)** in CSV format.

The input CSV file must contain the following columns:

| Column | Description | Unit |
|---------|-------------|------|
| `TRACK_ID` | Unique ID of each tracked particle | – |
| `POSITION_X` | X-coordinate | µm |
| `POSITION_Y` | Y-coordinate | µm |
| `FRAME` | Frame number | – |

---

# User-Defined Parameters

Several analysis parameters can be modified directly within each script.

Typical parameters include:

- Time interval between consecutive frames (`dt`)
- Frame rate (`frame_rate`)
- Frame interval used for velocity calculation (`frame_step`)
- Velocity threshold for Voronoi analysis
- Probe length range (`a`) with proper step size for Q₂ analysis
- Lag-time fitting range used to determine the anomalous diffusion exponent (α)
- Input and output file paths

Users should verify that these parameters are appropriate for their experimental conditions before running the scripts.

---

# Software Requirements

The scripts were developed using **Python 3.11**.

Required Python packages:

- numpy
- pandas
- matplotlib
- scipy
- shapely
- openpyxl

---

# Outputs

Depending on the script, the generated outputs include:

- High-resolution trajectory figures (.png)
- Velocity-colored trajectory plots
- Mean squared displacement (MSD) data
- Ensemble-averaged MSD
- Diffusion coefficients (D)
- Anomalous diffusion exponents (α)
- Displacement statistics
- Two-point overlap correlation functions Q₂(a,t)
- Voronoi seed points
- Voronoi cell areas
- Voronoi cell perimeters
- Voronoi neighbor counts
- Excel workbooks (.xlsx) containing all calculated quantities

---

# Notes

- These scripts were developed specifically for the analysis of SMT data presented in the accompanying manuscript.
- Users should modify only the input/output file paths and user-defined parameters appropriate for their datasets.
- All calculations are performed in two dimensions (2D).

---

# Citation

If you use these scripts in your research, please cite the associated manuscript:

> **...**

---
