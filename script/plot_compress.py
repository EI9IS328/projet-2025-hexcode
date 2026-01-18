#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt

out_dir = "script/img"
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv("results.csv")


methods = df["name"].unique()

plt.figure()

for method in methods:
    sub = df[df["name"] == method].sort_values("ex")
    plt.plot(
        sub["ex"],
        sub["writting_snapshots_time"],
        marker="o",
        label=method
    )

plt.xlabel("ex")
plt.ylabel("Writing time (ms)")
plt.title("Writting time relative to size")
plt.legend()
plt.grid(True)
plt.tight_layout()

out_path1 = os.path.join(out_dir, "writing_time.png")
plt.savefig(out_path1, dpi=300)
plt.close()





df_quant = df[df["name"] == "quant"].sort_values("ex")

plt.figure()

plt.plot(df_quant["ex"], df_quant["RMSE_moyen"], marker="o", label="avg_RMSE")
plt.plot(df_quant["ex"], df_quant["RMSE_max"], marker="s", label="max_RMSE")
plt.plot(df_quant["ex"], df_quant["Erreur_relative"], marker="^", label="relative_error")

plt.xlabel("ex")
plt.ylabel("RMSE / Relative error")
plt.title("RMSE and relative error compared to ex")
plt.legend()
plt.grid(True)
plt.tight_layout()

out_path2 = os.path.join(out_dir, "rmse_quant.png")
plt.savefig(out_path2)
plt.close()



df_quant = df[df["name"] == "quant"].sort_values("ex")

plt.figure()

plt.plot(df_quant["ex"], df_quant["RMSE_moyen"], marker="o", label="avg_RMSE")
plt.plot(df_quant["ex"], df_quant["RMSE_max"], marker="s", label="max_RMSE")

plt.xlabel("ex")
plt.ylabel("RMSE ")
plt.title("RMSE  compared to ex")
plt.legend()
plt.grid(True)
plt.tight_layout()

out_path2 = os.path.join(out_dir, "rmse_quant_without_rel.png")
plt.savefig(out_path2)
plt.close()



plt.figure()
for method in methods:
    sub = df[df["name"] == method].sort_values("ex")
    plt.plot(
        sub["ex"],
        sub["size_file_snapshots"],
        marker="o",
        label=method
    )

plt.xlabel("ex")
plt.ylabel("Size of snapshots (bytes)")
plt.title("Snapshots size relative to size")
plt.legend()
plt.grid(True)
plt.tight_layout()

out_path3 = os.path.join(out_dir, "snapshots_size.png")
plt.savefig(out_path3)
plt.close()


plt.figure()

sismos_methods = ["origin", "rle"]

for method in sismos_methods:
    sub = df[df["name"] == method].sort_values("ex")
    plt.plot(
        sub["ex"],
        sub["size_file_sismos"],
        marker="o",
        label=method
    )

plt.xlabel("ex")
plt.ylabel("Size of sismos (bytes)")
plt.title("Sismos size relative to ex")
plt.legend()
plt.grid(True)
plt.tight_layout()

# Sauvegarde
out_path4 = os.path.join(out_dir, "sismos_size.png")
plt.savefig(out_path4)
plt.close()


