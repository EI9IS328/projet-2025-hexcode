#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import sys 
import math

file = sys.argv[1]
df = pd.read_csv(file,sep=" ")  
ids = df["Index"].unique()

n_ids = len(ids)

# Disposition automatique (carrée)
cols = math.ceil(math.sqrt(n_ids))
rows = math.ceil(n_ids / cols)

plt.figure(figsize=(10,10))

for i, sensor_id in enumerate(ids, start=1):
    data = df[df["Index"] == sensor_id]
    
    plt.subplot(rows, cols, i)
    plt.plot(data["Pressure"], data["TimeStep"])
    plt.gca().invert_yaxis()
    plt.title(f"ID {sensor_id}")
    plt.xlabel("TimeStep")
    plt.ylabel("Pressure")
    plt.grid(True)

plt.tight_layout()
plt.show()