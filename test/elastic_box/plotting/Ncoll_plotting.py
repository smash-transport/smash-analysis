import os  # operating system interface
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.special as sp
import sys
from collections import OrderedDict
from glob import glob
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../../python_scripts')
from txt_io import save_table
import smash_basic_scripts as sb
import argparse

# Preambula configuration for plotting
execfile(os.path.dirname(os.path.abspath(__file__))+'/../../../python_scripts/common_plotting.py')

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--SMASH_data", required = True,
                    help = "File that contains extracted SMASH data. ")
parser.add_argument("--output", required = True,
                    help = "Path to the output file that will be created. ")
parser.add_argument("--setup", required = True,
                    help = "Test setup that is run. ")
parser.add_argument("--comp_prev_version", action='store_true',
                    help = "Plot comparison to previous SMASH version.")
args = parser.parse_args()

input_txt_files  = glob(args.SMASH_data)  # for geometic and stochastic criterion
output_pic_file = args.output
xvar            = args.setup.split('scatrate_vs_')[1] # Variable to be on x axis

# axis labels
x_labels = {'V': "Box volume V, fm$^3$",
            'sigma': "Elastic cross-section $\sigma$, fm$^2$",
            'N': "Particle number, $N$",
            'Ntest': "Testparticle factor, $N_{test}$",
            'T': "temperature $T$, GeV",
            'dt': "time step size $dt$, fm/c",
            'Kn': "Knudsen number"}

y_label =  ('${N_{coll}^{SMASH}}\\left/{N_{coll}^{theory}} \\right.$')
y_label_defintion =  ('$\\frac{N_{coll}^{SMASH}}{N_{coll}^{theory}} \, = \,$'
                      '$\\frac{1}{N \, N_{test} \,/ 2} \,'
                       '\\frac{N_{coll}}{N_{ev} t_{tot}}\,'
                       '\\frac{V}{\sigma N \, \langle v \\rangle}$')


def plot_data(input_txt_file, plot_position, plot_color):

    plt.subplot(plot_position)

    coll_criterion_name = os.path.basename(input_txt_file)[9:-4]
    plt.annotate(coll_criterion_name + " criterion",
                xy=(0.03, 0.075), xycoords='axes fraction',weight='heavy' , color=plot_color, fontsize =40)

    # Get data from file
    data = np.genfromtxt(input_txt_file, names=('Ncoll', 'Nevents',
                   't_run', 'V', 'sigma', 'N', 'Ntest', 'T', 'dt'))

    with open(input_txt_file, 'r') as f:
        smash_version = f.readlines()[-1].strip('# \n')

    # Sort by x axis variable
    data = data[data[xvar].argsort()]

    x = data[xvar]
    y = data['Ncoll']

    # Make a text label about the properties of the used box
    s=[]
    s.append('Elastic Box$:$')
    if (xvar != 'V'):     s.append("$V$ = %.1f fm$^3$"      % data['V'][0])
    if (xvar != 'sigma'): s.append("$\sigma$ = %.1f fm$^2$" % data['sigma'][0])
    if (xvar != 'N'):     s.append("$N$ = %i"               % data['N'][0])
    if (xvar != 'Ntest'): s.append("$N_{test}$ = %i"        % data['Ntest'][0])
    if (xvar != 'T'):     s.append("$T$ = %.3f GeV"         % data['T'][0])
    if (xvar != 'dt'):    s.append("$dt$ = %s fm/c"         % data['dt'][0])
    s.append("$t_{tot}$ = %.1f fm/c"    % data['t_run'][0])
    s.append("$N_{ev}$ = %i"            % data['Nevents'][0])
    box_label = '\n'.join(s)

    if plot_position == 211:  # only print title and input box once
        plt.annotate(box_label, xy=(1.02, 0.97), ha="left", va="top", xycoords='axes fraction', fontsize=30)
        plt.title('only $\pi^0$, only elastic collisions', fontsize=30, y=1.02)
    if plot_position == 212:
        plt.annotate(y_label_defintion, xy=(0.62, 0.125), xycoords='axes fraction', fontsize =40)


    # Number of collisions is expected to be equal to this norm (for <v> = c)
    norm = data['Nevents'] * data['t_run'] * (data['sigma'] * 0.5 * data['N'] * data['N'] * data['Ntest'] / data['V'])

    # Average relative velocity factor, arXiv:1311.4494, matters only at m/T < 0.7
    a = 0.135/data['T']  # m_pi0/T
    v_rel = 4./a * sp.kn(3, 2.0*a) / np.power(sp.kn(2, a), 2)
    norm *= v_rel

    y = y / norm
    y_error = np.sqrt(data['Ncoll']) / norm
    plt.errorbar(x, y, yerr=y_error, fmt='o', capsize=10,
                label=smash_version, markersize = 15,
                zorder = 2, markeredgecolor= plot_color, color=plot_color)

    if args.comp_prev_version:
        import comp_to_prev_version as cpv
        # plot reference data from previous SMASH version
        cpv.plot_previous_results('elastic_box', args.setup, '-' + coll_criterion_name + '.txt')

    plt.xlim(0.0, 1.05 * x.max())
    plt.ylim(ymin=0.4, ymax=max(1.5, 1.05 * y.max()))
    if (xvar == 'dt'):
        plt.xscale('log')
        plt.xlim(1.e-4, 2.0 * x.max())
        plt.gca().tick_params(pad=10)

    # store plotted data
    save_table(
        OrderedDict([('x', x), ('y', y), ("y_error", y_error)]),
        '{}.txt'.format(args.setup + "-" + coll_criterion_name),
        smash_version,
    )

    plt.legend(loc = 'upper right')
    plt.axhline(1, linewidth=3, linestyle='--', color='black', zorder = 0)
    plt.ylabel(y_label, fontsize=50)
    if plot_position == 212: plt.xlabel(x_labels[xvar])


plot_data(input_txt_files[0], 211, "midnightblue")
plot_data(input_txt_files[1], 212, "maroon")

plt.figtext(0.8, 0.95, "SMASH analysis: %s" % \
             (sb.analysis_version_string()), \
             color = "gray", fontsize = 10)
plt.tight_layout()

# Save picture and/or show it
plt.savefig(output_pic_file, bbox_inches = "tight", pad_inches=0.5)
#plt.show()
