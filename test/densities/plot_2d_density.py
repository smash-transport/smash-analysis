#!/usr/bin/env python3
import sys
from mpl_toolkits import mplot3d
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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
fig, axs = plt.subplots(nrows=nrows, ncols=ncols,figsize=(8,5))
fig.subplots_adjust(right=0.8, hspace=0.0, wspace=.05)
norm_mean = np.mean(Grid[0])
norm_min = np.min(Grid[0]) - .015 * norm_mean
norm_max = np.max(Grid[0]) + .015 * norm_mean
norm = mpl.colors.Normalize(vmin=norm_min, vmax=norm_max)

i = 0
for ax in axs.flat:
    tstep = times[i]
    rho_grid = Grid[i].reshape(dim_grid, dim_grid)
    ax.set_title('t = {} fm'.format(tstep), fontsize=15)
    im = ax.imshow(rho_grid, cmap='GnBu_r', interpolation='nearest', norm=norm, vmin = norm_min, vmax = norm_max)
    ax.set_xlabel('x', fontsize=18)
    if i > 0: ax.set_yticklabels([])
    else: ax.set_ylabel('y', fontsize=18)
    i = i + 1

plt.figtext(0.75, 0.8, " SMASH code:      %s\n SMASH analysis: %s\n" % \
            (smash_version, smash_analysis_version), \
            color = "gray", fontsize = 5)
cbar_ax = fig.add_axes([0.815, 0.236, 0.03, 0.526])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label(r'$\mathbf{\rho}$ [1/fm$^3$]',fontsize=16)
plt.savefig(outfile_name, bbox_inches='tight')
plt.close()
