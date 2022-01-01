# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import codecs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../python_scripts')
from collections import OrderedDict
from pdgs_from_config import pdg_to_name, pdg_to_name_init
import smash_basic_scripts as sb

parser = argparse.ArgumentParser()
parser.add_argument("--SMASH_data", required = True,
                    help = "File that contains extracted SMASH data. ")
parser.add_argument("--pdg1", required = True, type = int,
                    help = "PDG code of particle one.")
parser.add_argument("--pdg2", required = True, type = int,
                    help = "PDG code of particle two.")
parser.add_argument("--plab", required = True, type = float,
                    help = "Collision energy p_lab.")
parser.add_argument("--setup", required = True, type = str,
                    help = "Collision setup.")
parser.add_argument("--comp_prev_version", required = False, action='store_true',
                    help = "Plot comparison to previous SMASH version.")
args = parser.parse_args()


input_file = args.SMASH_data
pdg1 = args.pdg1
pdg2 = args.pdg2
p_lab = args.plab
setup = args.setup

base = os.path.splitext(input_file)[0]
output_file = base + ".pdf"

# Preambula configuration for plotting
exec(compile(open(os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts/common_plotting.py', "rb").read(), os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts/common_plotting.py', 'exec'))
colours = ["blue","red","green","darkmagenta","orange","deepskyblue", "magenta", "chartreuse", "cyan","limegreen"]
channels = ['total','N+N','N+N*','N*+Δ','N+Δ','N+Δ*','Δ+Δ','Δ+Δ*']
colour_coding = OrderedDict(list(zip(channels,colours)))

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

### (1) read particles and determine kinematics for plot
particle1 = pdg_to_name(int(pdg1),os.path.dirname(os.path.realpath(input_file)) + '/data/1/config.yaml')
particle2 = pdg_to_name(int(pdg2),os.path.dirname(os.path.realpath(input_file)) + '/data/1/config.yaml')
m1 = m2 = 0.938  # nucleon
sqrts = np.sqrt(m1**2 + m2**2 + 2 * m2 * np.sqrt(m1**2 + p_lab**2))

### (2) read reaction and SMASH version from header/last line
smash_version = codecs.open(input_file, 'r').readlines()[-1].split()[1]
with codecs.open(input_file, 'r', encoding='utf-8') as f:
    colnames = f.readline().rstrip().split()
ncols = len(colnames)

### (3) read data and divide by sin(theta)
contents = np.loadtxt(input_file, skiprows = 1, unpack=True)
theta = contents[0]

for i in range(0,len(theta)):
  s = np.sin(theta[i])
  if s>0.:
    for j in range(1, ncols):
      contents[j][i] = contents[j][i]/s

### (4) Plot different curves
for i in range(1, ncols):
  plt.plot(theta*180/np.pi, contents[i], label=str(colnames[i]), color=colour_coding.get(colnames[i]), zorder=2)


### (4a) plot data from previous SMASH version
if args.comp_prev_version:
    import comp_to_prev_version as cpv
    processes = ['total','N+N','N+N*','N*+Δ','N+Δ','N+Δ*','Δ+Δ','Δ+Δ*']
    cpv.plot_previous_results('angular_distributions', setup, '/theta.dat', color_list = colour_coding, process_list = processes)


### (5) set up axes, labels, etc
ax.set_xlim([0,180])
ax.set_xticks(list(range(0, 180, 10)), minor=True)
plt.title(particle1 + "+" + particle2 + " @ $\sqrt{s}=" + str(round(sqrts,2)) + "$ GeV")
plt.xlabel("$\\theta_{cm}$")
plt.ylabel("$d\sigma/dcos(\\theta_{cm})$ [mb]")
plt.yscale('log')
plt.figtext(0.77, 0.925, "SMASH analysis: %s" % \
             (sb.analysis_version_string()), \
             color = "gray", fontsize = 10)
plt.tight_layout()
plt.legend(title=smash_version, loc="best",fontsize = 20, ncol=2)

plt.savefig(output_file)
