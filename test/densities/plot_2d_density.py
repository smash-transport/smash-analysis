#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
import numpy as np
import math

fname = sys.argv[1]
outfile_name = sys.argv[2]

# Read information of header
with open(fname, 'r') as f:
    _, versions, \
    _, BoxLength, \
    _, dim_grid, \
    _, times = \
    f.readline(), f.readline(), \
    f.readline(), float(f.readline()), \
    f.readline(), int(f.readline()), \
    f.readline(), f.readline()

smash_version = versions.split()[0]
smash_analysis_version = versions.split()[1]
times = [float(i) for i in times.split()]
CellVol = (BoxLength / dim_grid) * (BoxLength / dim_grid) * BoxLength

data = np.loadtxt(fname, unpack=True, skiprows=8)

Grid = []
for num in range(len(times)):
    Grid.append(data[num] / CellVol)

ncols = 2 # fixed
nrows = int(math.ceil(float(len(times))/ncols))
fig, axs = plt.subplots(nrows=nrows, ncols=ncols,figsize=(8,6))
axs = axs.ravel()
fig.subplots_adjust(hspace=.1, wspace=.25)
norm_mean = np.mean(Grid[0])
norm_min = np.min(Grid[0]) - .015 * norm_mean
norm_max = np.max(Grid[0]) + .015 * norm_mean
norm = mpl.colors.Normalize(vmin=norm_min, vmax=norm_max)

for i, tstep in enumerate(times):
    rho_grid = Grid[i].reshape(dim_grid, dim_grid)
    axs[i].set_title('t = {} fm'.format(tstep), fontsize=15)
    im = axs[i].imshow(rho_grid, cmap='viridis', norm=norm)
    axs[i].set_xlabel('x', fontsize=18)
    axs[i].set_ylabel('y', fontsize=18)

plt.figtext(0.4, 0.8, " SMASH code:      %s\n SMASH analysis: %s\n" % \
            (smash_version, smash_analysis_version), \
            color = "gray", fontsize = 10)
cbar = fig.colorbar(im, ax=axs.ravel(), fraction=0.022, pad=0.04, norm=norm)
#cbar.set_label(r'$\langle N\rangle$',fontsize=18)
cbar.set_label(r'$\rho\;\rm [1 / fm^3]$',fontsize=16)
#plt.tight_layout()
plt.savefig(outfile_name)
