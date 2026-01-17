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

    plt.xlabel("Nombre d'éléments")
    plt.ylabel("Temps total (s)")
    plt.title("Temps total de la simulation in-situ en fonction du nombre d'éléments")
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
    SAVE_INTERVAL = 100
    df_filt = df[df["save_interval"] == SAVE_INTERVAL]

    plt.figure(figsize=(9, 6))

    plt.plot(df_filt["num_elements"]**3, df_filt["size_simsos_file"], marker="o", label="simsos")
    plt.plot(df_filt["num_elements"]**3, df_filt["size_slices_file"], marker="o", label="slices")
    plt.plot(df_filt["num_elements"]**3, df_filt["size_ppm_slices_file"], marker="o", label="ppm slices")

    plt.xlabel("Nombre d'éléments")
    plt.ylabel("Taille des fichiers (bytes)")
    plt.title(f"Taille des fichiers en fonction du nombre d'éléments\n(save_interval = {SAVE_INTERVAL})")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig("file_sizes.png")
    plt.show()

if __name__ == "__main__":
    main()
