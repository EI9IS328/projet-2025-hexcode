#!/usr/bin/env python3

import subprocess
import os
from glob import glob
import csv
import pandas as pd


def read_measures(measures_path):
    df = pd.read_csv(measures_path, sep=r'\s+')
    row = df.iloc[0]  

    return (
        float(row["kernel_time"]),
        float(row["size_file_snapshots"]),
        float(row["writting_snapshots_time"]),
        float(row.get("size_file_sismos", 0.0)),  
    )


def read_stat_compress(stat_path):
    df = pd.read_csv(stat_path, sep=r'\s+', header=None)
    row = df.iloc[0]  

    rmse_mean = float(row[0])
    rmse_max = float(row[1])
    rel_err = float(row[2])

    return rmse_mean, rmse_max, rel_err



def append_global_csv(
    csv_path,
    ex, ey, ez,
    name,
    measures_path,
    stat_path=None
):
    kernel_time, size_snap, write_time, size_sismos = read_measures(measures_path)

    if stat_path is not None and os.path.exists(stat_path):
        rmse_mean, rmse_max, rel_err = read_stat_compress(stat_path)
    else:
        rmse_mean = rmse_max = rel_err = ""

    with open(csv_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            ex, ey, ez, name,
            kernel_time,
            size_snap,
            write_time,
            size_sismos,
            rmse_mean,
            rmse_max,
            rel_err
        ])


def run_decompress(tab_path):
    for path in tab_path:
        subprocess.run(["python3", "./script/decompress.py", path])


def run_rmse(original_dirs, compressed_dirs):
    for orig, comp in zip(original_dirs, compressed_dirs):
        orig_snap = os.path.join(orig, "snapshots")
        comp_snap = os.path.join(comp, "snapshots")

        cmd = [
            "python3",
            "./script/rmse_quant.py",
            orig_snap,
            comp_snap
        ]

        subprocess.run(cmd, check=True)


def run_semproxy(val_range, extra_flags=None):
    if extra_flags is None:
        extra_flags = []

    created_dirs = []

    for val in val_range:
        cmd = [
            "./build/bin/semproxy",
            f"--ex={val}",
            f"--ey={val}",
            f"--ez={val}",
            "--mesh", "cartesian",
            "--save-snapshots",
            "--save-interval", "100"
        ] + extra_flags

        subprocess.run(cmd, check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        data_dirs = sorted(glob("data/data_*"))
        if data_dirs:
            latest_dir = data_dirs[-1]
            print(f"Dossier créé : {latest_dir}")
            created_dirs.append(latest_dir)

    return created_dirs


def main():
    val_range = range(10, 101, 10)

    original_dirs = run_semproxy(val_range, extra_flags=["--save-sismos"])
    compressed_quant_dirs = run_semproxy(val_range, extra_flags=["--quantify"])
    compressed_rle_dirs = run_semproxy(val_range, extra_flags=["--rle", "--save-sismos"])
    compressed_rle_quant_dirs = run_semproxy(val_range, extra_flags=["--rle", "--quantify"])
    
    run_decompress(compressed_quant_dirs)
    run_decompress(compressed_rle_dirs)
    run_decompress(compressed_rle_quant_dirs)

    run_rmse(original_dirs, compressed_quant_dirs)

    global_csv = "results.csv"

    with open(global_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "ex", "ey", "ez", "name",
            "kernel_time",
            "output_time",
            "traitement_time",
            "size_file_snapshots",
            "size_file_slices",
            "size_file_sismos",
            "writting_snapshots_time",
            "RMSE_moyen",
            "RMSE_max",
            "Erreur_relative"
        ])

    for val, orig, quant, rle, rle_quant in zip(
        val_range,
        original_dirs,
        compressed_quant_dirs,
        compressed_rle_dirs,
        compressed_rle_quant_dirs
    ):
        append_global_csv(
            csv_path=global_csv,
            ex=val, ey=val, ez=val,
            name="origin",
            measures_path=os.path.join(orig, "measure.csv"),
            stat_path=None
        )

        append_global_csv(
            csv_path=global_csv,
            ex=val, ey=val, ez=val,
            name="quant",
            measures_path=os.path.join(quant, "measure.csv"),
            stat_path=os.path.join(quant, "stat_compress.csv")
        )

        append_global_csv(
            csv_path=global_csv,
            ex=val, ey=val, ez=val,
            name="rle",
            measures_path=os.path.join(rle, "measure.csv"),
            stat_path=os.path.join(rle, "stat_compress.csv")
        )

        append_global_csv(
            csv_path=global_csv,
            ex=val, ey=val, ez=val,
            name="rle_quant",
            measures_path=os.path.join(rle_quant, "measure.csv"),
            stat_path=os.path.join(rle_quant, "stat_compress.csv")
        )


if __name__ == "__main__":
    main()
