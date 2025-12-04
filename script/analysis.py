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
    ap.add_argument('-d', '--data-folder', type=pathlib.Path, help="Folder containing the data")
    ap.add_argument('-o', '--output', type=str, default='analysis_results.csv',
                    help="Output CSV file name (default: analysis_results.csv)")
    args = ap.parse_args()

    mode = args.mode
    data_folder:pathlib.Path = args.data_folder
    outp = args.output_file

    snapshot_data = pd.DataFrame()
    if not data_folder.is_dir():
            print("Data folder must be a directory")
            sys.exit(1)
    snapshot_files = [f for f in (data_folder / 'snapshots').glob("snapshot_*") if f.is_file()]
    # sismo_files = [f for f in (data_folder / 'sismos').glob("sismo_*") if f.is_file()]

    if snapshot_files == []:
        print("No CSV files found in the specified directory")
        sys.exit(1)
    snapshot_data_frames = []
    for file in snapshot_files:
        df = pd.read_csv(file, sep=" ")
        snapshot_data_frames.append(df)

    snapshot_data = pd.concat(snapshot_data_frames, ignore_index=True)

    with open(outp, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        writer.writerow(["timestep", "min_pressure", "max_pressure", "mean_pressure", "median_pressure", "std_pressure"])

    if mode == 'global' or mode == 'all':
        with open(outp, 'a') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerow(["global",
                             snapshot_data['pressure'].min(),
                            snapshot_data['pressure'].max(),
                            snapshot_data['pressure'].mean(),
                            snapshot_data['pressure'].median(),
                            snapshot_data['pressure'].std()])

    if mode == 'step' or mode == 'all':
        timesteps = snapshot_data['timestep'].unique()
        for t in timesteps:
            timestep_data:pd.DataFrame = snapshot_data[snapshot_data['timestep'] == t]
            with open(outp, 'a') as f:
                writer = csv.writer(f, delimiter=' ')
                writer.writerow([t,
                                 timestep_data['pressure'].min(),
                                 timestep_data['pressure'].max(),
                                 timestep_data['pressure'].mean(),
                                 timestep_data['pressure'].median(),
                                 timestep_data['pressure'].std()])

if __name__ == "__main__":
    main()
