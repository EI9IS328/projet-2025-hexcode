#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
import pathlib
import csv

def load_dataframes_from_folder(folder:pathlib.Path, pattern:str):
    files = [f for f in folder.glob(pattern) if f.is_file()]
    if files == []:
        print(f"No files matching pattern {pattern} found in the specified directory")
        sys.exit(1)
    data_frames = []
    for file in files:
        df = pd.read_csv(file, sep=" ")
        data_frames.append(df)
    if data_frames:
        return pd.concat(data_frames, ignore_index=True)
    else:
        return pd.DataFrame()
    
def write_timestep_analysis(data:pd.DataFrame, outp):
    timesteps = data['timestep'].unique()
    for timestep in timesteps:
        timestep_data:pd.DataFrame = data[data['timestep'] == timestep]
        with open(outp, 'a') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerow([timestep,
                            data['pressure'].min(),
                            data['pressure'].max(),
                            data['pressure'].mean(),
                            data['pressure'].median(),
                            data['pressure'].std()])

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

    if not data_folder.is_dir():
            print("Data folder must be a directory")
            sys.exit(1)

    # Load data
    snapshot_data = pd.DataFrame()
    sismo_data = pd.DataFrame()
    if mode == 'global' or mode == 'step' or mode == 'all':
        snapshot_data = load_dataframes_from_folder(data_folder / 'snapshots', "snapshot_*")
    if mode == 'sismo' or mode == 'all':
        sismo_data = load_dataframes_from_folder(data_folder / 'sismos', "sismo_*")
    
    # Write headers
    with open(outp, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        if mode == 'all':
            writer.writerow(["timestep", "min_pressure", "max_pressure", "mean_pressure", "median_pressure", "std_pressure",
                         "sismo_min_pressure", "sismo_max_pressure", "sismo_mean_pressure", "sismo_median_pressure", "sismo_std_pressure"])
        else:
            writer.writerow(["timestep", "min_pressure", "max_pressure", "mean_pressure", "median_pressure", "std_pressure"])

    # Write global analysis
    if mode == 'global' or mode == 'all':
        with open(outp, 'a') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerow(["global",
                             snapshot_data['pressure'].min(),
                            snapshot_data['pressure'].max(),
                            snapshot_data['pressure'].mean(),
                            snapshot_data['pressure'].median(),
                            snapshot_data['pressure'].std()])

    # Write timestep analysis
    if mode == 'sismo' or mode == 'step' or mode == 'all':
        if mode == 'sismo':
            write_timestep_analysis(sismo_data, outp)
        elif mode == 'step':
            write_timestep_analysis(snapshot_data, outp)
        elif mode == 'all':
            timesteps = sismo_data['timestep'].unique()
            for timestep in timesteps:
                sismo_timestep_data:pd.DataFrame = sismo_data[sismo_data['timestep'] == timestep]
                snapshot_timestep_data:pd.DataFrame = snapshot_data[snapshot_data['timestep'] == timestep]
                with open(outp, 'a') as f:
                    writer = csv.writer(f, delimiter=' ')
                    writer.writerow([timestep,
                                    snapshot_timestep_data['pressure'].min(),
                                    snapshot_timestep_data['pressure'].max(),
                                    snapshot_timestep_data['pressure'].mean(),
                                    snapshot_timestep_data['pressure'].median(),
                                    snapshot_timestep_data['pressure'].std(),
                                    sismo_timestep_data['pressure'].min(),
                                    sismo_timestep_data['pressure'].max(),
                                    sismo_timestep_data['pressure'].mean(),
                                    sismo_timestep_data['pressure'].median(),
                                    sismo_timestep_data['pressure'].std()])

if __name__ == "__main__":
    main()
