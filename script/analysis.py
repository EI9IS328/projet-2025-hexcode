#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
import pathlib
import csv

# TODO : use pathlib for paths and csv module for csv writing

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-m', '--mode', type=str,
                    help="Range of analysis: 'global', 'step', 'sismo' or 'all'")
    ap.add_argument('-d', '--data-folder', type=str, help="Folder containing the data")
    ap.add_argument('-o', '--output', type=str, default='analysis_results.csv',
                    help="Output CSV file name (default: analysis_results.csv)")
    args = ap.parse_args()

    mode = args.mode
    data_folder = args.data_folder
    outp = args.output_file

    data = pd.DataFrame()

    if mode == 'global' or mode == 'all':
        import os
        if not os.path.isdir(data_folder):
            print("Data folder must be a directory")
            sys.exit(1)

        files = [f for f in os.listdir(data_folder) if f.endswith('.csv')]
        if files == []:
            print("No CSV files found in the specified directory")
            sys.exit(1)
        data_frames = []
        for file in files:
            df = pd.read_csv(os.path.join(data_folder, file), sep=" ")
            data_frames.append(df)

        data = pd.concat(data_frames, ignore_index=True)

    elif mode == 'step':
        data = pd.read_csv(data_folder, sep=" ")

    with open(outp, 'w') as f:
        f.write("timestep min_pressure max_pressure mean_pressure median_pressure std_pressure\n")

    if mode == 'global' or mode == 'all':
        with open(outp, 'a') as f:
            f.write(f"global {data['pressure'].min()} {data['pressure'].max()} "
                    f"{data['pressure'].mean()} {data['pressure'].median()} "
                    f"{data['pressure'].std()}\n")

    if mode == 'step' or mode == 'all':
        timesteps = data['timestep'].unique()
        for t in timesteps:
            timestep_data:pd.DataFrame = data[data['timestep'] == t]
            with open(outp, 'a') as f:
                f.write(f"{t} {timestep_data['pressure'].min()} {timestep_data['pressure'].max()} "
                        f"{timestep_data['pressure'].mean()} {timestep_data['pressure'].median()} "
                        f"{timestep_data['pressure'].std()}\n")    

if __name__ == "__main__":
    main()
