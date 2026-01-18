#!/usr/bin/env python3
import csv
import pathlib
import subprocess
import sys
import pandas as pd
import os

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "data"

NUM_ELEMENTS = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
SAVE_INTERVALS = [50, 100, 200]

# Save configurations: label -> CLI option
SAVE_CONFIGS = {
    "snapshots": ["--save-snapshots"],
    "sismos": ["--save-sismos"],
    "slices": ["--save-slices"],
    "ppm_slices": ["--save-ppm-slices"],
}


def find_latest_data(path):
    p = pathlib.Path(path)
    if not p.exists():
        return None
    folders = sorted(
        p.glob("data_*"),
        key=lambda x: x.stat().st_mtime,
        reverse=True
    )
    return folders[0] if folders else None


def parse_measure_file(path):
    data = pd.read_csv(path, sep=" ")
    return data.iloc[0]


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 benchmark_in-situ.py [OUTPUT_CSV]", file=sys.stderr)
        sys.exit(1)

    outp = sys.argv[1]
    sep = " "

    # Header
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

    for sizes in NUM_ELEMENTS:
        for interval in SAVE_INTERVALS:
            for save_type, save_flags in SAVE_CONFIGS.items():

                sem_proxy_cmd = [
                    str(PROJECT_ROOT / "build" / "bin" / "semproxy"),
                    "--ex", str(sizes),
                    "--ey", str(sizes),
                    "--ez", str(sizes),
                    "--in-situ",
                    "--save-interval", str(interval),
                    *save_flags
                ]

                try:
                    subprocess.run(
                        sem_proxy_cmd,
                        check=True,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL
                    )
                except subprocess.CalledProcessError as e:
                    print("semproxy failed:", e.returncode, file=sys.stderr)
                    sys.exit(1)

                data_path = find_latest_data(DEFAULT_DATA_DIR)

                if not data_path.exists() or data_path is None:
                    print("Data not found:", data_path, file=sys.stderr)
                    sys.exit(2)
                    

                measure_file = data_path / "measure.csv"
                if not measure_file.exists():
                    print("measure.csv not found in:", data_path, file=sys.stderr)
                    sys.exit(3)

                m = parse_measure_file(measure_file)

                # Raw times (µs)
                kernel_us = m["kernel_time"]
                output_us = m["output_time"]
                processing_us = m["traitement_time"]

                # Aggregated times (s)
                sem_seconds = (kernel_us + output_us) / 1e6
                analysis_seconds = processing_us / 1e6
                total_seconds = sem_seconds + analysis_seconds

                # File sizes (bytes)
                size_file_snapshots = m.get("size_file_snapshots", 0)
                size_file_slices = m.get("size_file_slices", 0)
                size_file_sismos = m.get("size_file_sismos", 0)
                size_file_ppm_slices = m.get("size_file_ppm_slices", 0)
                size_file_analysis = m.get("size_file_analysis", 0)
                # Write result
                with open(outp, "a") as f:
                    writer = csv.writer(f, delimiter=sep)
                    writer.writerow([
                        "in-situ",
                        sizes,
                        interval,
                        save_type,
                        kernel_us,
                        output_us,
                        processing_us,
                        sem_seconds,
                        analysis_seconds,
                        total_seconds,
                        size_file_snapshots,
                        size_file_slices,
                        size_file_sismos,
                        size_file_ppm_slices,
                        size_file_analysis
                    ])


if __name__ == "__main__":
    main()
