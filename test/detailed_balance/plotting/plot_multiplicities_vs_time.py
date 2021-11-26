import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../../python_scripts')
import smash_basic_scripts as sb

input_file = sys.argv[1]
output_file = sys.argv[2]
tstart = float(sys.argv[3])
config_file = sys.argv[4]

# Preambula configuration for plotting
exec(compile(open(os.path.dirname(os.path.abspath(__file__))+'/../../../python_scripts/common_plotting.py', "rb").read(), os.path.dirname(os.path.abspath(__file__))+'/../../../python_scripts/common_plotting.py', 'exec'))

with open(input_file, 'r') as f:
    f.readline()
    smash_version = f.readline()
    f.readline()
    total_events = int(f.readline())
    f.readline()
    pdg_list = f.readline().rstrip().split()

contents = np.loadtxt(input_file, skiprows=6)
time = contents[0]
n_pdgs = len(pdg_list)

plt.xlabel("t [fm/c]", fontsize = 36)
plt.ylabel("number of particles", fontsize = 36)

plt.title("Box: chemical equilibration")
smash_code_version, analysis_version = smash_version.split()
plt.figtext(0.8, 0.94, " %d events\n SMASH code:      %s\n SMASH analysis: %s" % \
             (total_events, smash_code_version, analysis_version), \
             color = "gray", fontsize = 10)

for i in np.arange(n_pdgs):
    label = sb.pdg_to_name(int(pdg_list[i]), config_file)
    plt.plot(time, contents[2*i + 1], label=label, lw=3)
    x = time[time > tstart]
    y = contents[2*i + 1][time > tstart]
    A = np.polyfit(x, y, 0)
    plt.axhline(y=A[0], xmin=tstart/x.max(),
                linewidth=1, color='k', linestyle='-')
plt.axvline(x=tstart, linewidth=3, color='k', linestyle='--')

plt.legend(ncol = 3, fontsize = 30)

plt.savefig(output_file)
