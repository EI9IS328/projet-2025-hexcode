#!/usr/bin/env python3
import csv
import pathlib
import subprocess
import sys
import time
import pandas as pd
import os

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "data"

NUM_ELEMENTS = [i for i in range(10, 81, 10)]
SAVE_INTERVALS = [50, 100, 200]
SAVE_TYPES = ["snapshots","sismos"]

def find_latest_data(path):
    p = pathlib.Path(path)
    if not p.exists():
        return None
    folder = sorted(p.glob("data_*"), key=lambda x: x.stat().st_mtime, reverse=True)
    return folder[0] if folder else None

def main():
    if sys.argv[1] is None:
        print("Usage: python3 benchmark.py [OUTPUT_CSV]", file=sys.stderr)
        sys.exit(1)

    outp = sys.argv[1]
    sep = " "
    is_empty = not os.path.exists(outp) or os.path.getsize(outp) == 0

    with open(outp, "a", newline="") as f:
        writer = csv.writer(f, delimiter=sep)

        if is_empty:
            writer.writerow([
                "approche",
                "num_elements",
                "save_interval",
                "save_type",
                "kernel_us",
                "output_us",
                "processing_us",
                "sem_seconds",
                "analysis_seconds",
                "total_seconds",
                "size_file_snapshots",
                "size_file_slices",
                "size_file_sismos",
                "size_file_ppm_slices",
                "size_file_analysis"
            ])

    for size in NUM_ELEMENTS:
        for save_type in SAVE_TYPES:
            for interval in SAVE_INTERVALS:
                # sismos does not depend on save_intervals
                if save_type!="sismos" or (save_type=="sismos" and interval==SAVE_INTERVALS[0]):
                    sem_proxy_cmd = [str(PROJECT_ROOT / "build" / "bin" / "semproxy"),
                                    "--ex", str(size),
                                    "--ey", str(size),
                                    "--ez", str(size),
                                    "--mesh", "cartesian",
                                    "--save-"+save_type,
                                    "--save-interval", str(interval)]

                    # Run semproxy
                    try:
                        subprocess.run(sem_proxy_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    except subprocess.CalledProcessError as e:
                        print("semproxy failed with return code", e.returncode, file=sys.stderr)
                        sys.exit(1)

                    data_path = find_latest_data(DEFAULT_DATA_DIR)

                    if not data_path.exists() or data_path is None:
                        print("Data not found:", data_path, file=sys.stderr)
                        sys.exit(2)
                    
                    if not (data_path / "measure.csv").exists():
                        print("Measure file not found in:", data_path, file=sys.stderr)
                        sys.exit(3)

                    measures = pd.read_csv(data_path / "measure.csv", sep=' ')

                    kernel_us = measures["kernel_time"][0]
                    output_us = measures["output_time"][0]
                    processing_us = measures["traitement_time"][0]

                    sem_seconds = (kernel_us + output_us + processing_us) / 1e6

                    # Run analysis.py
                    analysis_cmd = [sys.executable,
                                    str(PROJECT_ROOT / "script" / "analysis.py"),
                                    "-m", save_type,
                                    "-d", str(data_path),
                                    "-o", str(data_path / "analysis_results.csv")]
                    t0 = time.perf_counter()
                    try:
                        subprocess.run(analysis_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    except subprocess.CalledProcessError as e:
                        print("analysis.py failed with return code", e.returncode, file=sys.stderr)
                        sys.exit(4)
                    t1 = time.perf_counter()
                    analysis_seconds = t1 - t0

                    total_seconds = sem_seconds + analysis_seconds

                    size_analysis = (data_path / "analysis_results.csv").stat().st_size

                    # Write to CSV
                    with open(outp, "a") as f:
                        writer = csv.writer(f, delimiter=sep)
                        writer.writerow(["ad-hoc",size, interval, save_type,
                            kernel_us, output_us, processing_us, sem_seconds, analysis_seconds, total_seconds,
                            measures["size_file_snapshots"][0], measures["size_file_slices"][0], measures["size_file_sismos"][0],measures["size_file_ppm_slices"][0], size_analysis])


if __name__ == "__main__":
    main()
