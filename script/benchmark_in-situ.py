#!/usr/bin/env python3
import csv
import pathlib
import subprocess
import sys
import time
import pandas as pd

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = "data"

NUM_ELEMENTS = [10,20,30,40,50,60,70,80,90,100]
SAVE_INTERVALS = [50, 100, 200]

def find_latest_data(path):
    p = pathlib.Path(path)
    if not p.exists():
        return None
    folder = sorted(p.glob("data_*"), key=lambda x: x.stat().st_mtime, reverse=True)
    return folder[0] if folder else None

def parse_measure_file(path):
    measures = {}
    with open(path, "r") as f:
        data = pd.read_csv(f,sep=" ")
        for key in data.columns:
            measures[key] = data[key].tolist() 
        pass
    return measures

def main():
    if len (sys.argv) < 2:
        print("Usage: python3 benchmark_in-situ.py [OUTPUT_CSV]", file=sys.stderr)
        sys.exit(1)

    outp = sys.argv[1]
    sep = " "
    with open(outp, "w") as f:
        writer = csv.writer(f, delimiter=sep)
        writer.writerow(["num_elements", "save_interval", "sem_seconds", "analysis_seconds", "total_seconds","size_simsos_file","size_snapshots_file","size_slices_file","size_ppm_slices_file"])

    for sizes in NUM_ELEMENTS:
        for interval in SAVE_INTERVALS:
            sem_proxy_cmd = [
                             str(PROJECT_ROOT / "build" / "bin" / "semproxy"),
                             "--ex", str(sizes),
                             "--ey", str(sizes),
                             "--ez", str(sizes),
                             "--save-snapshots",
                             "--save-sismos",
                             "--save-slices",
                             "--save-ppm-slices",
                             "--in-situ",
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

            measures = parse_measure_file(data_path / "measure.csv")

            kernel_ms = measures.get("kernel_time")[0]
            output_ms = measures.get("output_time")[0]
            traitement_ms = measures.get("traitement_time")[0]

            sem_seconds = (kernel_ms + output_ms) / 1e6
            analysis_seconds = traitement_ms / 1e6
            total_seconds = sem_seconds + analysis_seconds
            
            size_simsos_file = measures.get("size_file_sismos")[0]
            size_snapshots_file = measures.get("size_file_snapshots")[0]
            size_slices_file = measures.get("size_file_slices")[0]
            size_ppm_slices_file = measures.get("size_file_ppm_slices")[0]

            # Write to CSV
            with open(outp, "a") as f:
                writer = csv.writer(f, delimiter=sep)
                writer.writerow([sizes, interval, sem_seconds, analysis_seconds, total_seconds, size_simsos_file, size_snapshots_file, size_slices_file, size_ppm_slices_file])


if __name__ == "__main__":
    main()