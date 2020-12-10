import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append(
    os.path.dirname(
        os.path.abspath(__file__)) +
    '/../../python_scripts')
from read_binary_output import *
import smash_basic_scripts as sb
import glob
import yaml

# Preambula configuration for plotting
execfile(os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts/common_plotting.py')

DataOutfile = sys.argv[1]
PlotName = sys.argv[2]
try:
  comp_prev_version = sys.argv[3]
  comp_prev_version = True
except IndexError:
  comp_prev_version = False

# Read information of header
with open(DataOutfile, 'r') as f:
    _, _, versions = \
    f.readline(), f.readline(), f.readline()
smash_version = versions.split()[1]
smash_analysis_version = versions.split()[2]

TimeSteps, dens_arr = np.loadtxt(DataOutfile, unpack=True, skiprows=5)

plt.plot(TimeSteps, dens_arr, label=smash_version, color = 'darkred')

# old version ?
if comp_prev_version:
    import comp_to_prev_version as cpv
    processes = []
    cpv.plot_previous_results('densities', DataOutfile.split('/')[-2], '/' + DataOutfile.split('/')[-1])

plt.title(r'Density central cell', loc='left', fontsize = 30)
plt.figtext(0.8, 0.95, "SMASH analysis: %s" % \
             (sb.analysis_version_string()), \
             color = "gray", fontsize = 10)
plt.xlabel(r't [fm]')
plt.ylabel(r'$\rho_\mathsf{B}/\rho_0$')
plt.xlim(0,40)
plt.ylim(0,5)
plt.tight_layout()
plt.legend()
plt.savefig(PlotName, bbox_inches = "tight")
plt.close()
