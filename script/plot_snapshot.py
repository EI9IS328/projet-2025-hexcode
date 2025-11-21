#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns




file = sys.argv[1] #file name with the time step to analyse

data = pd.read_csv(file,sep=" ")
tstep =  data.at[1,'timestep']
yz = data[(data['x'] == 100) & (data['timestep'] == tstep)]
xy = data[(data['z'] == 100) & (data['timestep'] == tstep)]
xz = data[(data['y'] == 100) & (data['timestep'] == tstep)]

tmp = yz.pivot(index = 'y',columns='z',values='pressure')
plt.subplot(221)
sns.heatmap(tmp,vmin = yz['pressure'].min(),vmax=yz['pressure'].max(),cmap='seismic')
tmp = xy.pivot(index = 'x',columns='y',values='pressure')
plt.subplot(222)
sns.heatmap(tmp,vmin = yz['pressure'].min(),vmax=yz['pressure'].max(),cmap='seismic')
tmp = xz.pivot(index = 'x',columns='z',values='pressure')
plt.subplot(223)
sns.heatmap(tmp,vmin = yz['pressure'].min(),vmax=yz['pressure'].max(),cmap='seismic')


plt.show()
