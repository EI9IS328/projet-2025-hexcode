#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import matplotlib.ticker as mticker



file = sys.argv[1] 

data = pd.read_csv(file,sep=" ")
tstep =   str(data.at[1,'timestep'])

yz = data[(data['plane'] == 'yz') ]
xy = data[(data['plane'] == 'xy') ]
xz = data[(data['plane'] == 'xz') ]

tmp = xy.pivot(index = 'i',columns='j',values='pressure')
plt.subplot(221)
sns.heatmap(tmp,vmin = xy['pressure'].min(),vmax=xy['pressure'].max(),cmap='seismic')
plt.title("XY slice")
ax = plt.gca()
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
tmp = xz.pivot(index = 'i',columns='j',values='pressure')
plt.subplot(222)
sns.heatmap(tmp,vmin = xz['pressure'].min(),vmax=xz['pressure'].max(),cmap='seismic')
plt.title("XZ slice")
ax = plt.gca()
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
tmp = yz.pivot(index = 'i',columns='j',values='pressure')
plt.subplot(223)
sns.heatmap(tmp,vmin = yz['pressure'].min(),vmax=yz['pressure'].max(),cmap='seismic')
plt.title("YZ slice")
ax = plt.gca()
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))

plt.suptitle("Pressure Field 2D Slices - " + tstep)

plt.subplots_adjust(wspace=0.5, hspace=0.6)
plt.show()