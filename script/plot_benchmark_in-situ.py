#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_results.py merged.csv")
        sys.exit(1)

    csv_file = sys.argv[1]

    # Lecture du CSV
    df = pd.read_csv(csv_file, sep=",")

    SAVE_INTERVAL = 100

    # Filtrage in-situ + save_interval
    df_in_situ = df[
        (df["approche"] == "in-situ") &
        (df["save_interval"] == SAVE_INTERVAL)
    ]

    # ==================================================
    # STACKPLOT 1 — SISMOS (in-situ)
    # ==================================================
    df_sismos = df_in_situ[df_in_situ["save_type"] == "sismos"].sort_values("num_elements")

    x = df_sismos["num_elements"]**3
    kernel = df_sismos["kernel_us"] / 1e6
    analysis = df_sismos["analysis_seconds"]

    plt.figure(figsize=(9, 6))

    plt.stackplot(
        x,
        kernel,
        analysis,
        labels=["Kernel", "Analysis"],
        alpha=0.8
    )

    plt.xlabel("Nombre d'éléments³")
    plt.ylabel("Temps (s)")
    plt.title(f"Décomposition du temps — Sismos")
    plt.legend(loc="upper left")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("stackplot_sismos_in_situ.png")
    plt.show()

    # ==================================================
    # STACKPLOT 2 — SNAPSHOTS (in-situ)
    # ==================================================
    df_snap = df_in_situ[df_in_situ["save_type"] == "snapshots"].sort_values("num_elements")

    x = df_snap["num_elements"]**3
    kernel = df_snap["kernel_us"] / 1e6
    analysis = df_snap["analysis_seconds"]

    plt.figure(figsize=(9, 6))

    plt.stackplot(
        x,
        kernel,
        analysis,
        labels=["Kernel", "Analysis"],
        alpha=0.8
    )

    plt.xlabel("Nombre d'éléments³")
    plt.ylabel("Temps (s)")
    plt.title(f"Décomposition du temps — Snapshots")
    plt.legend(loc="upper left")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("stackplot_snapshots_in_situ.png")
    plt.show()
    
    SAVE_INTERVAL = 100
    df = df[df["save_interval"] == SAVE_INTERVAL]

    # ===============================
    # -------- GRAPHE 1 : SISMOS
    # ===============================

    df_sismos = df[df["save_type"] == "sismos"]

    df_sismos_in = df_sismos[df_sismos["approche"] == "in-situ"].sort_values("num_elements")
    df_sismos_ad = df_sismos[df_sismos["approche"] == "ad-hoc"].sort_values("num_elements")

    # Temps total (secondes)
    df_sismos_in["total_time"] = (
        df_sismos_in["analysis_seconds"] +
        df_sismos_in["kernel_us"] / 1e6
    )

    df_sismos_ad["total_time"] = (
        df_sismos_ad["sem_seconds"] +
        df_sismos_ad["analysis_seconds"]
    )

    plt.figure(figsize=(9, 6))

    plt.plot(
        df_sismos_in["num_elements"]**3,
        df_sismos_in["total_time"],
        marker="o",
        linewidth=2,
        label="sismos in-situ"
    )

    plt.plot(
        df_sismos_ad["num_elements"]**3,
        df_sismos_ad["total_time"],
        marker="x",
        linestyle="--",
        linewidth=2,
        label="sismos ad-hoc"
    )

    plt.xlabel("Nombre d'éléments³")
    plt.ylabel("Temps total (s)")
    plt.title(f"Temps total – Sismos")

    plt.ylim(
        0,
        max(df_sismos_in["total_time"].max(),
            df_sismos_ad["total_time"].max()) * 1.1
    )

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # ===============================
    # -------- GRAPHE 2 : SNAPSHOTS
    # ===============================

    df_snap = df[df["save_type"] == "snapshots"]

    df_snap_in = df_snap[df_snap["approche"] == "in-situ"].sort_values("num_elements")
    df_snap_ad = df_snap[df_snap["approche"] == "ad-hoc"].sort_values("num_elements")

    # Temps total (secondes)
    df_snap_in["total_time"] = (
        df_snap_in["analysis_seconds"] +
        df_snap_in["kernel_us"] / 1e6
    )

    df_snap_ad["total_time"] = (
        df_snap_ad["sem_seconds"] +
        df_snap_ad["analysis_seconds"]
    )

    plt.figure(figsize=(9, 6))

    plt.plot(
        df_snap_in["num_elements"]**3,
        df_snap_in["total_time"],
        marker="o",
        linewidth=2,
        label="snapshots in-situ"
    )

    plt.plot(
        df_snap_ad["num_elements"]**3,
        df_snap_ad["total_time"],
        marker="x",
        linestyle="--",
        linewidth=2,
        label="snapshots ad-hoc"
    )

    plt.xlabel("Nombre d'éléments³")
    plt.ylabel("Temps total (s)")
    plt.title(f"Temps total – Snapshots")

    plt.ylim(
        0,
        max(df_snap_in["total_time"].max(),
            df_snap_ad["total_time"].max()) * 1.1
    )

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
