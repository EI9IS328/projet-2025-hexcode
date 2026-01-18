#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_ad_hoc_perf.py result.csv")
        sys.exit(1)

    csv_file = sys.argv[1]

    df = pd.read_csv(csv_file, sep=" ")
    snap_df = df[df["save_type"] == "snapshots"]
    sismo_df = df[df["save_type"] == "sismos"]
    sismo_df = sismo_df[sismo_df["save_interval"] == 100]

    # snapshot write time
    plt.figure(figsize=(8, 6))

    for interval in sorted(snap_df["save_interval"].unique()):
        sub = snap_df[snap_df["save_interval"] == interval]
        plt.plot(
            sub["num_elements"],
            sub["output_us"]/1e6,
            marker="o",
            label=f"save_interval={interval}"
        )

    plt.xlabel("Nombre d'éléments (par dimension)")
    plt.ylabel("Temps total (s)")
    plt.title("Temps d'écritures des snapshots en fonction du nombre\nd'éléments et de l'intervalle de sauvegarde des snapshots")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("snapshot_write_time.png")

    # snapshot write size
    plt.figure(figsize=(8, 6))

    for interval in sorted(snap_df["save_interval"].unique()):
        sub = snap_df[snap_df["save_interval"] == interval]
        plt.plot(
            sub["num_elements"],
            sub["size_file_snapshots"],
            marker="o",
            label=f"save_interval={interval}"
        )

    plt.xlabel("Nombre d'éléments (par dimension)")
    plt.ylabel("Taille totales des fichiers (en octets)")
    plt.title("Taille combinée de l'ensemble des snapshots écrites durant la simulation")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("snapshot_write_size.png")


    # snap total
    plt.figure(figsize=(8, 6))

    sub = snap_df[snap_df["save_interval"] == 100]

    plt.stackplot(
        sub["num_elements"],
        sub["kernel_us"]/1e6,
        sub["output_us"]/1e6,
        sub["analysis_seconds"],
        labels=["Temps du kernel", "Temps d'écritures", "Temps d'analyses"],
        alpha=0.8
    )

    plt.xlabel("Nombre d'éléments (par dimension)")
    plt.ylabel("Temps total (s)")
    plt.title("Différence entre les temps d'écritures et les temps d'analyses des\nsnapshots et l'exécution du kernel en fonction du nombre d'éléments")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("total_snapshots_time.png")


    # sismo total
    plt.figure(figsize=(8, 6))

    plt.stackplot(
        sismo_df["num_elements"],
        sismo_df["kernel_us"]/1e6,
        sismo_df["output_us"]/1e6,
        sismo_df["analysis_seconds"],
        labels=["Temps kernel", "Temps écritures", "Temps analyses"],
        alpha=0.8
    )

    plt.xlabel("Nombre d'éléments (par dimension)")
    plt.ylabel("Temps total (s)")
    plt.title("Différence entre les temps d'écritures et les temps d'analyses des\nsismogrammes et l'exécution du kernel en fonction du nombre d'éléments")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("total_sismos_time.png")


if __name__ == "__main__":
    main()