#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import sys
import pathlib
import subprocess
import time
import csv 

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1] 
SEM_PROXY = str(PROJECT_ROOT / "build" / "bin" / "semproxy")
DATA_DIR = PROJECT_ROOT / "data"
DECOMPRESS_SCRIPT = PROJECT_ROOT / "script" / "decompress.py"

NUM_ELEMENTS = [10, 15, 20, 25, 30]
SAVE_INTERVALS = [50, 100, 200]

def find_latest_data(path):
    p = pathlib.Path(path)
    if not p.exists():
        return None
    folders = sorted(p.glob("data_*"),key=lambda x: x.stat().st_mtime,reverse=True)
    return folders[:2] if folders else None

def main():
    outp = sys.argv[1]
    results = []
    for sizes in NUM_ELEMENTS:
        for interval in SAVE_INTERVALS:
            sem_proxy_cmd_original = [
                             str(SEM_PROXY),
                             "--ex", str(sizes),
                             "--ey", str(sizes),
                             "--ez", str(sizes),
                             "--mesh", "cartesian",
                             "--save-snapshots",
                             "--save-interval", str(interval)]
            sem_proxy_cmd_compress = [
                             str(SEM_PROXY),
                             "--ex", str(sizes),
                             "--ey", str(sizes),
                             "--ez", str(sizes),
                             "--mesh", "cartesian",
                             "--save-snapshots",
                             "--save-interval", str(interval),
                             "--compress"]
            t0 = time.perf_counter()
            subprocess.run(sem_proxy_cmd_original, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            t1 = time.perf_counter()
            time_save_original = t1 -t0
            
            t0 = time.perf_counter()

            subprocess.run(sem_proxy_cmd_compress, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            t1 = time.perf_counter()
            time_save_compress = t1 - t0
            
            dir = find_latest_data(DATA_DIR)
            original_dir = dir[1]
            recon_dir = dir[0]
            rmse_list = []
            rel_max_list = []
            t0 = time.perf_counter()
            subprocess.run(
                    [sys.executable,str(DECOMPRESS_SCRIPT), str(recon_dir / "snapshots")],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    cwd=str(PROJECT_ROOT)
                )
            t1 = time.perf_counter()
            decompress_time = t1 - t0
            
            for filename in sorted(os.listdir(original_dir / "snapshots")):
                orig_path = os.path.join(original_dir / "snapshots", filename)
                recon_path = os.path.join(recon_dir / "snapshots", filename)
                df_orig = pd.read_csv(orig_path, sep=r"\s+", engine="python")
                df_rec  = pd.read_csv(recon_path, sep=r"\s+", engine="python")
                df_orig.columns = df_orig.columns.str.strip()
                df_rec.columns  = df_rec.columns.str.strip()
                p_orig = df_orig["pressure"].to_numpy(dtype=float)
                p_rec  = df_rec["pressure"].to_numpy(dtype=float)
                error = p_rec - p_orig
                rmse = np.sqrt(np.mean(error**2))
                rmse_list.append(rmse)
                mask = np.abs(p_orig) > 1e-12
                if np.any(mask):
                    rel_max = np.max(np.abs(error[mask] / p_orig[mask]))
                else:
                    rel_max = 0.0  
                rel_max_list.append(rel_max)
                print(f"{filename} → RMSE = {rmse:.6e}, RelMax = {rel_max:.6e}")
            
            
            results.append({
                "time_process_original" : time_save_original,
                "time_process_compress" : time_save_compress,
                "time_process_decompress" : decompress_time,
                "RMSE avg" : np.mean(rmse_list),
                "RMSE max" : np.max(rmse_list),
                "Max relative error" : np.max(rel_max_list)
            })
            df = pd.DataFrame(results)
            df.to_csv(outp)

if __name__ == "__main__":
    main()




          