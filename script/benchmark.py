#!/usr/bin/env python3
import csv
import pathlib
import subprocess
import sys
import time

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1] / "projet-2025-hexcode"
DEFAULT_DATA_DIR = PROJECT_ROOT / "data"

NUM_ELEMENTS = [5, 10, 15, 20, 25, 30]
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
        # TODO
        pass
    return measures

def main():
    if sys.argv[1] is None:
        print("Usage: python3 benchmark.py [OUTPUT_CSV]", file=sys.stderr)
        sys.exit(1)

    outp = sys.argv[1]
    sep = " "
    with open(outp, "w") as f:
        writer = csv.writer(f, delimiter=sep)
        writer.writerow(["num_elements", "save_interval", "sem_seconds", "analysis_seconds", "total_seconds"])

    for sizes in NUM_ELEMENTS:
        for interval in SAVE_INTERVALS:
            sem_proxy_cmd = [sys.executable,
                             str(PROJECT_ROOT / "build" / "bin" / "semproxy"),
                             "--ex", str(sizes),
                             "--ey", str(sizes),
                             "--ez", str(sizes),
                             "--mesh", "cartesian",
                             "-s",
                             "--save-sismos",
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
            
            if not (data_path / "mesure").exists():
                print("Mesure file not found in:", data_path, file=sys.stderr)
                sys.exit(3)

            measures = parse_measure_file(data_path / "")

            kernel_ms = measures.get("kernel")[0]
            output_ms = measures.get("output")[0]

            sem_seconds = (kernel_ms + output_ms) / 1e6

            # Run analysis.py
            analysis_cmd = [sys.executable,
                            str(PROJECT_ROOT / "script" / "analysis.py"),
                            "-m", "all",
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

            # Write to CSV
            with open(outp, "a") as f:
                writer = csv.writer(f, delimiter=sep)
                writer.writerow([sizes, interval, sem_seconds, analysis_seconds, total_seconds])


if __name__ == "__main__":
    main()