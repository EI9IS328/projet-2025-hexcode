#!/usr/bin/env python3
import pandas as pd
import sys

if sys.argv[1] is None or sys.argv[2] is None:
    print("Usage: python3 analysis.py <mode> <data_folder or filename> [output_file]")
    sys.exit(1)

mode = sys.argv[1] #range of analysis
data_folder = sys.argv[2] #folder name or filename
out = sys.argv[3] if len(sys.argv) > 3 else "analysis_results.csv"

if mode not in ['global', 'step', "all"]:
    print("Mode must be 'global', 'step' or 'all'")
    sys.exit(1)

data = pd.DataFrame()

if mode == 'global' or mode == 'all':
    import os
    if not os.path.isdir(data_folder):
        print("Data folder must be a directory")
        sys.exit(1)
    
    files = [f for f in os.listdir(data_folder) if f.endswith('.csv')]

    data_frames = []
    for file in files:
        df = pd.read_csv(os.path.join(data_folder, file), sep=" ")
        data_frames.append(df)

    data = pd.concat(data_frames, ignore_index=True)

elif mode == 'step':
    data = pd.read_csv(data_folder, sep=" ")

with open(out, 'w') as f:
    f.write("timestep min_pressure max_pressure mean_pressure median_pressure std_pressure\n")

if mode == 'global' or mode == 'all':
    with open(out, 'a') as f:
        f.write(f"global {data['pressure'].min()} {data['pressure'].max()} "
                f"{data['pressure'].mean()} {data['pressure'].median()} "
                f"{data['pressure'].std()}\n")

if mode == 'step' or mode == 'all':
    timesteps = data['timestep'].unique()
    for t in timesteps:
        timestep_data:pd.DataFrame = data[data['timestep'] == t]
        with open(out, 'a') as f:
            f.write(f"{t} {timestep_data['pressure'].min()} {timestep_data['pressure'].max()} "
                    f"{timestep_data['pressure'].mean()} {timestep_data['pressure'].median()} "
                    f"{timestep_data['pressure'].std()}\n")    
