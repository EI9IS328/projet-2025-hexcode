#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import sys 
import math

file = sys.argv[1]
df = pd.read_csv(file,sep=" ")  
ids = df["Index"].unique()

n_ids = len(ids)

cols = math.ceil(math.sqrt(n_ids))
rows = math.ceil(n_ids / cols)

plt.figure(figsize=(10,10))

for i, sensor_id in enumerate(ids, start=1):
    data = df[df["Index"] == sensor_id]
    coordX = int(df.loc[df["Index"] == sensor_id , "X"].iloc[0])
    coordY = int(df.loc[df["Index"] == sensor_id , "Y"].iloc[0])
    coordZ = int(df.loc[df["Index"] == sensor_id , "Z"].iloc[0])

    plt.subplot(rows, cols, i)
    plt.plot(data["Pressure"], data["TimeStep"])
    plt.gca().invert_yaxis()
    plt.title(f"({coordX},{coordY},{coordZ})Recv {sensor_id}")
    plt.xlabel("TimeStep")
    plt.ylabel("Pressure")
    plt.grid(True)

plt.suptitle("Pressure evolution for different recever")
plt.tight_layout()
plt.show()