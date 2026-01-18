#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import sys 
import math
import os

file = sys.argv[1]
df = pd.read_csv(file,sep=" ")  
ids = df["index"].unique()
out_dir = "script/img"
os.makedirs(out_dir, exist_ok=True)

n_ids = len(ids)

cols = math.ceil(math.sqrt(n_ids))
rows = math.ceil(n_ids / cols)

plt.figure(figsize=(10,10))

for i, sensor_id in enumerate(ids, start=1):
    data = df[df["index"] == sensor_id]
    coordX = int(df.loc[df["index"] == sensor_id , "x"].iloc[0])
    coordY = int(df.loc[df["index"] == sensor_id , "y"].iloc[0])
    coordZ = int(df.loc[df["index"] == sensor_id , "z"].iloc[0])

    plt.subplot(rows, cols, i)
    plt.plot(data["pressure"], data["timestep"])
    plt.gca().invert_yaxis()
    plt.title(f"({coordX},{coordY},{coordZ})Recv {sensor_id}")
    plt.xlabel("timestep")
    plt.ylabel("pressure")
    plt.grid(True)

plt.suptitle("Pressure evolution for different recever")
out_path4 = os.path.join(out_dir, "sismo.png")
plt.savefig(out_path4)