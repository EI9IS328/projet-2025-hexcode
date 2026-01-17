#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_results.py result.csv")
        sys.exit(1)

    csv_file = sys.argv[1]

    # Lecture du CSV (séparateur = espace)
    df = pd.read_csv(csv_file, sep=" ")

    # =========================
    # 1. Temps d'exécution
    # =========================
    plt.figure(figsize=(8, 6))

    for interval in sorted(df["save_interval"].unique()):
        sub = df[df["save_interval"] == interval]
        plt.plot(
            sub["num_elements"]**3,
            sub["total_seconds"],
            marker="o",
            label=f"save_interval={interval}"
        )

    plt.xlabel("Num elements")
    plt.ylabel("Total time (s)")
    plt.title("Total execution time vs problem size")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("total_time.png")
    plt.show()

    # =========================
    # 2. Détail SEM vs Analyse
    # =========================
    plt.figure(figsize=(8, 6))

    interval = df["save_interval"].unique()[0]
    sub = df[df["save_interval"] == interval].sort_values("num_elements")

    x = sub["num_elements"]**3
    sem = sub["sem_seconds"]
    analysis = sub["analysis_seconds"]

    plt.stackplot(
        x,
        sem,
        analysis,
        labels=["SEM time", "Analysis time"],
        alpha=0.8
    )

    plt.xlabel("Nombre d'éléments")
    plt.ylabel("Temps (s)")
    plt.title(f"Total du temps de la simulation in-situ (save_interval={interval})")
    plt.legend(loc="upper left")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("time_breakdown_stacked.png")
    plt.show()

    # =========================
    # 3. Taille des fichiers
    # =========================
    plt.figure(figsize=(8, 6))

    plt.plot(df["num_elements"]**3, df["size_simsos_file"], marker="o", label="Sismos")
    plt.plot(df["num_elements"]**3, df["size_snapshots_file"], marker="o", label="Snapshots")
    plt.plot(df["num_elements"]**3, df["size_slices_file"], marker="o", label="Slices")
    plt.plot(df["num_elements"]**3, df["size_ppm_slices_file"], marker="o", label="PPM Slices")

    plt.xlabel("Nombre d'éléments")
    plt.ylabel("Taille des fichiers")
    plt.title("Taille des fichiers de sortie vs taille du problème")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("file_sizes.png")
    plt.show()


if __name__ == "__main__":
    main()
