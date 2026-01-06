#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import sys

original_dir = sys.argv[1]
recon_dir = sys.argv[2]
rmse_list = []
rel_max_list = []

for filename in sorted(os.listdir(original_dir)):
    if not filename.endswith(".csv"):
        continue

    orig_path = os.path.join(original_dir, filename)
    recon_path = os.path.join(recon_dir, filename)

    if not os.path.exists(recon_path):
        print(f"Fichier manquant : {filename}")
        continue

    df_orig = pd.read_csv(orig_path, sep=r"\s+", engine="python")
    df_rec  = pd.read_csv(recon_path, sep=r"\s+", engine="python")

    df_orig.columns = df_orig.columns.str.strip()
    df_rec.columns  = df_rec.columns.str.strip()

    if "pressure" not in df_orig.columns:
        raise RuntimeError(f"Colonne 'pressure' absente dans {filename}")

    p_orig = df_orig["pressure"].to_numpy(dtype=float)
    p_rec  = df_rec["pressure"].to_numpy(dtype=float)

    error = p_rec - p_orig

    rmse = np.sqrt(np.mean(error**2))
    rmse_list.append(rmse)

    mask = np.abs(p_orig) > 1e-12
    rel_max = np.max(np.abs(error[mask] / p_orig[mask]))
    rel_max_list.append(rel_max)

    print(f"{filename} → RMSE = {rmse:.6e}, RelMax = {rel_max:.6e}")

print("\n===== STATISTIQUES GLOBALES =====")
print(f"RMSE moyen          : {np.mean(rmse_list):.6e}")
print(f"RMSE max            : {np.max(rmse_list):.6e}")
print(f"Erreur relative max : {np.max(rel_max_list):.6e}")

