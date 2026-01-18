#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
import pathlib
import csv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-m', '--mode', type=str,
                    help="Range of analysis: 'snapshot' or 'sismo'", default='snapshot')
    ap.add_argument('-d', '--data-folder', type=pathlib.Path, help="Folder containing the data")
    ap.add_argument('-o', '--output', type=str, default='analysis_results.csv',
                    help="Output CSV file name (default: analysis_results.csv)")
    args = ap.parse_args()

    mode = args.mode
    data_folder:pathlib.Path = args.data_folder
    outp = args.output

    if not data_folder.is_dir():
            print("Data folder must be a directory")
            sys.exit(1)

    data_files = []
    if mode == 'snapshot':
        data_files = [f for f in (data_folder / 'snapshots').glob("snapshot_*") if f.is_file()]
    else:
        data_files = [f for f in (data_folder / 'sismos').glob("sismo") if f.is_file()]

    if data_files == []:
        print("No CSV files found in the specified directory")
        sys.exit(1)
    
    data_frames = []
    for file in data_files:
        df = pd.read_csv(file, sep=" ")
        data_frames.append(df)

    data:pd.DataFrame = pd.concat(data_frames, ignore_index=True)

    with open(outp, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        if mode == 'snapshot':
            writer.writerow(["timestep", "min", "max", "mean", "median", "std"])
        else:
            writer.writerow(["receiver_index", "min", "max", "mean", "median", "std"])

    if mode == 'snapshot':
        with open(outp, 'a') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerow(["global",
                             data['pressure'].min(),
                             data['pressure'].max(),
                             data['pressure'].mean(),
                             data['pressure'].median(),
                             data['pressure'].std()])
    index = []
    if mode == 'sismo':
        index = data['index'].unique()
    else:
        index = data['timestep'].unique()
    for i in index:
        row_data:pd.DataFrame = data[data['timestep'] == i] if mode == 'snapshot' else data[data['index'] == i]
        with open(outp, 'a') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerow([i,
                            row_data['pressure'].min(),
                            row_data['pressure'].max(),
                            row_data['pressure'].mean(),
                            row_data['pressure'].median(),
                            row_data['pressure'].std()])

if __name__ == "__main__":
    main()
