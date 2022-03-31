import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import argparse  # command line argument parser
from definitions import *  # IMPORT BINS, LINESTYLES AND CHANNELS
from glob import glob
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../python_scripts')
import smash_basic_scripts as sb

matplotlib.rcParams.update({'font.size': 16,
                            'axes.labelsize': 22,
                            'mathtext.default': 'rm',
                            'legend.fontsize': 15})

desc = """
        Plotting script for dilepton analysis.

        """

parser = argparse.ArgumentParser(description=desc)
parser.add_argument("--system", help="collision system that was used")
parser.add_argument("--energy", required = True,
                    help="kinetic energy (Ekin) per nucleon of target")
parser.add_argument("--data_dir", nargs='?', default="", required = False,
                    help="directory containing the experimental data")
parser.add_argument("--comp_prev_version", action = 'store_true',
                    help = "Plot comparison to previous SMASH version.")
args = parser.parse_args()

plot_with_data = False
if args.data_dir != "":
    plot_with_data = True
    print("Plotting with (HADES) data ...")

cross_sections_dict = {"pp1.25": 46.96,
                       "pp2.2": 42.16,
                       "pp3.5": 43.40,
                       "pNb3.5": 848.0}

def version():
    data = np.genfromtxt("other.version.dat", dtype='str')
    return data[1]

def normalization_AA():  # CC and ArKCl
    with open("other.avg_pion.dat") as avg:
        raw_norm = np.loadtxt(avg)
        norm_AA = raw_norm[1]  # use average of (n_piz+n_pim)/2
        print("Using avg. no. of pion =", norm_AA, "for normalization off AA spectra ...")
    return norm_AA

# create wider bins


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

def store_results(results_file, x, y, version):
    store_file = open(results_file, "w")
    store_file.write("m" + "\t" + "dN/dm" + "\n")
    for i in range(0,len(x)):
        store_file.write(str(x[i]) + "\t" + str(y[i]) + "\n")
    store_file.write("# " + str(version))

def plot(name, bin_factor, ch_list, style_dict, datafile="", cut_legend=""):

    n_ch = len(ch_list)

    # import hist_data
    with open("hist_raw_" + name + ".txt") as df_smash:
        data = np.loadtxt(df_smash, delimiter=' ', unpack=True)

    bin_centers = data[0, :]
    hist = data[1:n_ch + 1, :]

    # make dN/dx plot
    bin_width = bin_centers[1] - bin_centers[0]
    hist_dx = hist[:] / bin_width

    # renormalize for data comparison (currently only mass spectra compared with data)
    if datafile != "":
        # do cross section plot for pp, data in mub
        if args.system == "pp" or args.system == "pNb":
            hist_dx = hist_dx * \
                cross_sections_dict[args.system + args.energy] * 1000
        # spectra for CC is normalized with averaged number of pions
        if args.system == "CC" or args.system == "ArKCl":
            hist_dx = hist_dx / normalization_AA()

    # rebin
    bin_centers_new, hist_new = rebin(bin_centers, hist_dx, n_ch, bin_factor)

    # header for final hist. files
    head = []
    head.append(name)

    # plotting
    plt.plot(bin_centers_new, sum(hist_new),
             label="all", color='k', linewidth=3)
    for i in range(len(ch_list)):
        if (sum(hist_new[i]) > pow(10,-8)):  # any contributions to channel?
            plt.plot(bin_centers_new,
                     hist_new[i], style_dict["l_style"][i], label=ch_list[i], linewidth=2)
        head.append(ch_list[i])
    # save final hist. files
    np.savetxt("dN_d" + name +".txt", np.transpose(np.concatenate((bin_centers_new[None,:], hist_new))), header='    '.join(head))
    
    # plot data
    if datafile != "":
        if 'mass' in name:
            data_path = os.path.join(args.data_dir + 'm_inv_spectrum_dileptons/' + args.system ,datafile)
        elif 'pt' in name:
            data_path = os.path.join(args.data_dir + 'pT_spectrum/' + args.system + '/',datafile)
        elif 'y' in name:
            data_path = os.path.join(args.data_dir + 'y_spectrum/' + args.system + '/',datafile)
        else:
            print('No experimental data found.')

        with open(data_path) as df:
            data = np.loadtxt(df, unpack=True)

        x_data = data[0, :]
        y_data = data[1, :]
        y_data_err = data[2, :]

        plt.errorbar(x_data, y_data, yerr=y_data_err,
                     fmt='ro', ecolor='k', label="HADES", zorder = 3)

    # write and read old results, only for total dN/dm spectra
    if name == "mass":
        # store dN/dm spectra for future comparison to previous version
        store_results("dN_dm_tot.txt", bin_centers_new, sum(hist_new), version())

        if args.comp_prev_version:
            import comp_to_prev_version as cpv
            # plot reference data from previous SMASH version
            setup = str(args.system) + "_" + str(args.energy) + "_" + "filtered"
            cpv.plot_previous_results('dileptons', setup, '/dN_dm_tot.txt')


    leg_cols = 2
    if ("rho" in name) or ("omega" in name): leg_cols = 1

    # plot style
    leg = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=leg_cols, mode="expand",
                     borderaxespad=0, fancybox=True, title=args.system + "@" + args.energy + " GeV " + cut_legend)
    #plt.annotate(version(), xy=(0.02, 0.03), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), fontsize=12)
    plt.annotate('SMASH version: ' + version() + '\n' +
                'Analysis suite version: ' + sb.analysis_version_string(),
                xy=(0.4, 0.92), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"), fontsize=8)
    plt.xlim(style_dict["x_min"], style_dict["x_max"])
    plt.ylim(style_dict["y_min"], style_dict["y_max"])
    plt.xlabel(style_dict["xlab"])
    plt.ylabel(style_dict["ylab"])
    plt.yscale('log')
    plt.savefig("plot_" + name + ".pdf",
                bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.cla()


# PLOTS #

# rebinning factors (no rebinning for now)
mass_bf = 0
pt_bf = 0
rap_bf = 0

# mass spectra where data is available
if plot_with_data:
    if args.system == "pp":
        plot("mass",       mass_bf, ch_list_main, style_dict_mass_w_data_pp_pNb,
             datafile="ekin" + args.energy + ".exp")
    if args.system == "pNb":
        plot("mass",       mass_bf, ch_list_main, style_dict_mass_w_data_pp_pNb,
             datafile="ekin" + args.energy + ".exp")
        plot("mass_0_800", mass_bf, ch_list_main, style_dict_mass_w_data_pp_pNb,
             datafile="ekin" + args.energy + "_m0_m800.exp")
        plot("mass_800",   mass_bf, ch_list_main, style_dict_mass_w_data_pp_pNb,
             datafile="ekin" + args.energy + "_m800.exp")
    if args.system == "CC" or args.system == "ArKCl":
        plot("mass",       mass_bf, ch_list_main, style_dict_mass_w_data_CC_ArKCl,
             datafile="ekin" + args.energy + ".exp")
else:
    plot("mass", mass_bf, ch_list_main, style_dict_mass)


# pt spectra where data is available
if plot_with_data and args.system == "pp" and args.energy == "3.5":
        plot("pt_0_150",     pt_bf,   ch_list_main, style_dict_pt_w_data,
             datafile="ekin" + args.energy + "_m0_m150_dileptons.exp", cut_legend=", m < 150 MeV")
        plot("pt_150_470",   pt_bf,   ch_list_main, style_dict_pt_w_data,
             datafile="ekin" + args.energy + "_m150_m470_dileptons.exp", cut_legend=", 150 MeV < m < 470 MeV")
        plot("pt_470_700",   pt_bf,   ch_list_main, style_dict_pt_w_data,
             datafile="ekin" + args.energy + "_m470_m700_dileptons.exp", cut_legend=", 470 MeV < m < 700 MeV")
        plot("pt_700",       pt_bf,   ch_list_main, style_dict_pt_w_data,
             datafile="ekin" + args.energy + "_m700_dileptons.exp", cut_legend=", 700 MeV < m")
else:
    plot("pt_0_150",   pt_bf,   ch_list_main,  style_dict_pt, cut_legend=", m < 150 MeV")
    plot("pt_150_470", pt_bf,   ch_list_main,  style_dict_pt, cut_legend=", 150 MeV < m < 470 MeV")
    plot("pt_470_700", pt_bf,   ch_list_main,  style_dict_pt, cut_legend=", 470 MeV < m < 700 MeV")
    plot("pt_700",     pt_bf,   ch_list_main,  style_dict_pt, cut_legend=", 700 MeV < m")


# y spectra where data is available
if plot_with_data and args.system == "pp" and args.energy == "3.5":
        plot("y_0_150",      rap_bf,  ch_list_main,  style_dict_y_w_data,
             datafile="ekin" + args.energy + "_m0_m150_dileptons.exp", cut_legend=", m < 150 MeV")
        plot("y_150_470",    rap_bf,  ch_list_main,  style_dict_y_w_data,
             datafile="ekin" + args.energy + "_m150_m470_dileptons.exp", cut_legend=", 150 MeV < m < 470 MeV")
        plot("y_470_700",    rap_bf,  ch_list_main,  style_dict_y_w_data,
             datafile="ekin" + args.energy + "_m470_m700_dileptons.exp", cut_legend=", 470 MeV < m < 700 MeV")
        plot("y_700",        rap_bf,  ch_list_main,  style_dict_y_w_data,
             datafile="ekin" + args.energy + "_m700_dileptons.exp", cut_legend=", 700 MeV < m")
else:
    plot("y_0_150",   rap_bf,   ch_list_main,  style_dict_y, cut_legend=", m < 150 MeV")
    plot("y_150_470", rap_bf,   ch_list_main,  style_dict_y, cut_legend=", 150 MeV < m < 470 MeV")
    plot("y_470_700", rap_bf,   ch_list_main,  style_dict_y, cut_legend=", 470 MeV < m < 700 MeV")
    plot("y_700",     rap_bf,   ch_list_main,  style_dict_y, cut_legend=", 700 MeV < m")


# spectra where no data is available
plot("pt",             pt_bf,   ch_list_main,  style_dict_pt)
plot("y",       rap_bf,  ch_list_main,  style_dict_y)

# origin plots
plot("mass_rho",       mass_bf, ch_list_rho,   style_dict_mass_origin)
plot("mass_omega",     mass_bf, ch_list_omega, style_dict_mass_origin)

plot("pt_rho",         pt_bf,   ch_list_rho,   style_dict_pt_origin)
plot("pt_omega",       pt_bf,   ch_list_omega, style_dict_pt_origin)

plot("y_rho",   rap_bf,  ch_list_rho,   style_dict_y_origin)
plot("y_omega", rap_bf,  ch_list_omega, style_dict_y_origin)
