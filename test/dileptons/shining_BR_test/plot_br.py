import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import sys
import argparse  # command line argument parser
matplotlib.rcParams.update({'font.size': 16,
                            'axes.labelsize': 22,
                            'mathtext.default': 'rm',
                            'legend.fontsize': 15})

desc = """
        Modified plotting script from dilepton analysis
        for use of plotting the rho BR of pions vs
        dileptons scaled by the BR of the dilepton decay.

        """

# disabled error message, since there will be some bins
np.seterr(all='ignore')
# with no contribution, which would produce errors for log-scaling and ratios


parser = argparse.ArgumentParser(description=desc)
parser.add_argument("system", help="collision system that was used")
parser.add_argument(
    "energy", help="kinetic energy (Ekin) per nucleon of target")
args = parser.parse_args()


# channels to analyze
ch_list = [r'$\rho \rightarrow e^+e^-$',
           r'$\omega \rightarrow e^+e^-$',
           r'$\phi \rightarrow e^+e^-$',
           r'$\pi^0 \rightarrow \gamma e^+e^-$',
           r'$\eta \rightarrow \gamma e^+e^-$',
           r"$\eta' \rightarrow \gamma e^+e^-$",
           r'$\omega \rightarrow \pi^0 e^+e^-$',
           r'$\phi \rightarrow \pi^0 e^+e^-$',
           r'$\Delta^0\rightarrow n e^+e^-$',
           r'$\Delta^+ \rightarrow p e^+e^-$',
           r'unknown']
n_ch = len(ch_list)

colors = ['b', 'g', 'r', 'c', 'm', 'y']
nc = len(colors)
linestyles = ['-', '--', '-.', ':']


def rebin(x, y, ch, bin_factor):

    if bin_factor == 0:
        return x, y

    cut = len(x) % bin_factor
    if cut > 0:
        x = x[:-cut]
        y = y[:, :-cut]

    x_new_list = []
    y_new_list = [[] for i in range(ch)]

    for i in range(0, len(x), bin_factor):
        x_new_list.append(sum(x[i: i + bin_factor]) / bin_factor)
        for c in range(ch):
            y_new_list[c].append(sum(y[c][i: i + bin_factor]) / bin_factor)

    x_new = np.asarray(x_new_list)
    len(x_new)
    y_new = [[] for i in range(ch)]
    for c in range(ch):
        y_new[c] = np.asarray(y_new_list[c])
        len(y_new[c])

    return x_new, y_new


# the bin_factor increases the bin size (by default bin size is 1 MeV, which is usally to small)
def plots(spectrum, bin_factor):

    with open("hist_" + spectrum + "_11.txt") as datafile:
        data = np.loadtxt(datafile, delimiter=' ', unpack=True)

    x = data[0, :]
    y = data[1:n_ch + 1, :]

    x_ee, y_ee = rebin(x, y, n_ch, bin_factor)

    with open("hist_mass_211.txt") as datafile:
        data = np.loadtxt(datafile, delimiter=' ', unpack=True)

    x = data[0, :]
    y = data[1:n_ch + 1, :]

    x_pipi, y_pipi = rebin(x, y, n_ch, bin_factor)

    BR_pi = y_pipi[0] / (y_pipi[0] + y_ee[0])
    BR_e = y_ee[0] / (y_pipi[0] + y_ee[0])

    plt.plot(x_pipi, BR_pi, 'b-', label=r'$B.R.(\rho \rightarrow \pi^+\pi^-)$')
    plt.plot(x_ee, BR_e / (4.72e-5), 'r-',
             label=r'$B.R.(\rho \rightarrow e^+e^-)/4.72e-5$')

    # plot styling
    if spectrum == "mass":
        xaxis = r'$m\,[GeV]$'
        yaxis = r'$BR(m)$'
        plt.ylim(1E-3, 1E4)
        plt.xlim(0.0, 1.2)
        # leg = plt.legend(bbox_to_anchor=(0.7, 1.05), loc='upper left', ncol=1, fancybox=True)

    # good default legend (more specific legends commented out, can be used for nicer  manual plots)
    leg = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand",
                     borderaxespad=0, fancybox=True, title=args.system + "@" + args.energy + " GeV")

    plt.axvline(x=0.775, alpha=0.5, linestyle="--", color="k")  # rho pole mass

    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.yscale('log')
    plt.savefig("plot_" + spectrum + "_br.pdf",
                bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.cla()


# bin factors
mass_bf = 10

plots("mass", mass_bf)
