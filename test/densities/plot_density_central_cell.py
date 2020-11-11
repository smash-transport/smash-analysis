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

matplotlib.rcParams.update({'font.size': 16,
                            'axes.labelsize': 22,
                            'mathtext.default':'rm',
                            'legend.fontsize': 10,
                            'font.sans-serif':['Helvetica'],
                            'font.family':'sans-serif'})

DataOutfile = sys.argv[1]
PlotName = sys.argv[2]
try:
  comp_prev_version = sys.argv[3]
  comp_prev_version = True
except IndexError:
  comp_prev_version = False

if comp_prev_version:
    import comp_to_prev_version as cpv
    processes = []
    cpv.plot_previous_results('densities', setup, '/' + DataOutfile.split('/')[-1])

# Read information of header
with open(DataOutfile, 'r') as f:
    _, _, versions = \
    f.readline(), f.readline(), f.readline()
smash_version = versions.split()[1]
smash_analysis_version = versions.split()[2]

TimeSteps, dens_arr = np.loadtxt(DataOutfile, unpack=True, skiprows=5)

# old version ?
plt.plot(TimeSteps, dens_arr, linewidth=2, linestyle="-", label='AuAu')
# plt.text(0.3, 4.5, " SMASH code:      %s\n SMASH analysis: %s\n" % \
#             (smash_version, smash_analysis_version), \
#             color = "gray", fontsize = 7)
plt.title(r'Density central cell', loc='left')
plt.figtext(0.58, 0.9, " SMASH code:      %s\n SMASH analysis: %s" % \
             (smash_version, smash_analysis_version), \
             color = "gray", fontsize = 7)
plt.xlabel(r'$t\rm [fm]$')
plt.ylabel(r'$\rho_B/\rho_0$')
plt.xlim(0,40)
plt.ylim(0,5)
plt.tight_layout()
plt.legend()
plt.savefig(PlotName)
