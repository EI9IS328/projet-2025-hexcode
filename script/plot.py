#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys


file = sys.argv[1]
nRecv = int(sys.argv[2])

data = pd.read_csv(file,sep=" ")
data2 = data.loc[data["nb"]==nRecv,:]
plt.plot(data2["time"],data2["pression"])
plt.xlabel("Time(ms)")
plt.ylabel("Pression")
plt.title("Pression calculé en fonction du temps")
plt.savefig('zae.png', bbox_inches='tight')

