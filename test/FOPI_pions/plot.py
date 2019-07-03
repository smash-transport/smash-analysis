#!/usr/bin/python2
#
# Plots the pion multiplicity and ratio from SMASH versus the FOPI data.

import sys  # system-specific functions
import os  # operating system interface
from collections import OrderedDict
from glob import glob
import numpy as np
import matplotlib
import argparse
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../python_scripts')
from txt_io import save_table
from common_plotting import smash_style
import smash_basic_scripts as sb
import comp_to_prev_version as cpv


# Parse arguments from command line (experimental data if available)
parser = argparse.ArgumentParser()
parser.add_argument("--FOPI_multiplicity", required = False, default = '',
                    help = "File that contains FOPI multiplicity data. ")
parser.add_argument("--FOPI_ratio", required = False, default = '',
                    help = "File that contains FOPI ratio data. ")
parser.add_argument("--comp_prev_version", action='store_true',
                    help = "Plot comparison to previous SMASH version.")
args = parser.parse_args()

### (1) read SMASH data
names = ('e',
         'Npiplus', 'Npiplus_err',
         'Npizero', 'Npizero_err',
         'Npiminus','Npiminus_err',
         'n_event', 'n_test',
         'pt_piplus', 'pt_piplus_err',
         'pt_pizero', 'pt_pizero_err',
         'pt_piminus','pt_piminus_err')

data_file = os.getcwd() + "/pion_yields.txt"
smash_data = np.loadtxt(data_file, unpack=True)
smash = {}
for i, d in enumerate(smash_data):
    smash[names[i]] = d

with open(data_file, 'r') as f:
    version = f.readline().rstrip().split()[1]     # read SMASH version

### (2) plot multiplicity

plt.figure(figsize=(15, 30))
plt.figtext(0.79, 0.96, "SMASH analysis: %s" % \
             (sb.analysis_version_string()), \
             color = "gray", fontsize = 8)

plt.subplot(311)
plt.title("Au + Au (b < 2 fm)")

if args.FOPI_multiplicity != '':
    import comp_to_exp_data as ced
    ced.plot_FOPI_data(args.FOPI_multiplicity, 'multiplicity')

x = smash['e']
N = smash['n_test']
y = 1.5*(smash['Npiplus'] + smash['Npiminus'])
yerr = np.sqrt(y)/N
y /= N
plt.errorbar(x, y, yerr=yerr, label=version, color="darkred", linestyle = '-', zorder=3)

# store results
save_table(OrderedDict([('x', x), ('y', y), ('yerr', yerr)]), 'multiplicity.txt', version)

if args.comp_prev_version: # Compare to previous version
    import comp_to_prev_version as cpv
    cpv.plot_previous_results('FOPI_pions', 'multiplicity', '.txt')

smash_style.set(line_styles=False, update_legends=True)

plt.xlabel("$E_{kin}$ [AGeV]")
plt.xlim(0.3, 1.6)
plt.yscale('log')
plt.ylim(2, 100)
plt.ylabel(r"$M(\pi)=3/2\,(N_{\pi^+}+N_{\pi^-})$")
plt.legend(loc="lower right")

### (3) plot ratio

plt.subplot(312)
if args.FOPI_ratio != '':
    if 'ced' not in sys.modules:
        import comp_to_exp_data as ced
    ced.plot_FOPI_data(args.FOPI_ratio, 'ratio')

x = smash['e']
y = smash['Npiminus']/smash['Npiplus']
yerr = y * np.sqrt((smash['Npiminus_err'] / smash['Npiminus'])**2 +
                   (smash['Npiplus_err'] / smash['Npiplus'])**2)
plt.errorbar(x, y, yerr=yerr, label=version, linestyle = '-', color = "darkred",
             zorder = 3)

# store results
save_table(OrderedDict([('x', x), ('y', y), ('yerr', yerr)]), 'ratio.txt', version)

if args.comp_prev_version: # Compare to previous version
    # Compare to previous version
    cpv.plot_previous_results('FOPI_pions', '', 'ratio.txt')

plt.xlabel("$E_{kin}$ [AGeV]")
plt.xlim(0.3, 1.6)
plt.ylabel("$N_{\pi^-}/N_{\pi^+}$")
plt.legend(loc="best")

### (4) plot pt

plt.subplot(313)

plt.xlim(0.3, 1.6)
plt.xlabel('$E_{kin}$ [AGeV]')
plt.ylabel(r'$\left<p_T\right>$ [GeV]')
plt.errorbar(smash['e'], smash['pt_piplus'], yerr = smash['pt_piplus_err'],
            label='$\\pi^+$', color = "midnightblue" )
plt.errorbar(smash['e'], smash['pt_piminus'], yerr = smash['pt_piminus_err'],
            label='$\\pi^-$', color = "orange")
plt.legend(loc='best')

smash_style.set(line_styles=False, update_legends=True)
plt.tight_layout()

plt.savefig("FOPI_pions.pdf")
