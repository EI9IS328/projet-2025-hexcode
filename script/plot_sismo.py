#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import sys 

file = sys.argv[1]
df = pd.read_csv(file,sep=" ")  
ids = df["Index"].unique()

plt.figure(figsize=(14,6))

# Boucle pour tracer une courbe par id
for sensor_id in ids:
    data = df[df["Index"] == sensor_id]
    plt.plot(data["TimeStep"], data["Pressure"], label=f"ID {sensor_id}")

plt.title("Pression sismique par capteur")
plt.xlabel("Temps")
plt.ylabel("Pression")
plt.grid(True)
plt.legend(title="Capteurs")
plt.tight_layout()
plt.show()