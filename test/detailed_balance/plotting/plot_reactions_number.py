#/usr/bin/python
# coding=UTF-8

import sys
import os  # operating system interface
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../../python_scripts')
import smash_basic_scripts as sb
from glob import glob
from itertools import cycle
import matplotlib.gridspec as gridspec
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--SMASH_data", required = True,
                    help = "File that contains extracted SMASH data. ")
parser.add_argument("--output", required = True,
                    help = "Path to the output file that will be created. ")
parser.add_argument("--config", required = True,
                    help = "Path to the config.yaml file. ")
parser.add_argument("--setup", required = True,
                    help = "Test setup that is run. ")
parser.add_argument("--comp_prev_version", required = False, action='store_true',
                    help = "Plot comparison to previous SMASH version.")
args = parser.parse_args()

input_file = args.SMASH_data
output_file = args.output
config_file = args.config
setup = args.setup
results_file = os.path.dirname(args.SMASH_data)+'/Nreact_by_Nisopingroup.txt'

reload(sys)
sys.setdefaultencoding('utf-8')

def isospin_form(react):
    lhs = sorted([sb.pdg_to_name(x ,config_file).translate(None, '0+-⁺⁰⁻\^{}$').replace('p','n') for x in react[0]])
    rhs = sorted([sb.pdg_to_name(x, config_file).translate(None, '0+-⁺⁰⁻\^{}$').replace('p','n') for x in react[1]])
    return ' '.join(lhs) + '-' + ' '.join(rhs)

def clebsch_from_str(s):
    """For strings like 7/2 returns 3.5, for strings like '1' returns 1.0"""
    if ('/' in s):
       return float(s.split('/')[0])/float(s.split('/')[1])
    else:
       return float (s)

def get_descr(clebsch, part_in, part_out, arrow):
    """Return text description of a reaction in->out."""
    descr = u''
    for i in xrange(part_in.shape[0]):
        descr += sb.pdg_to_name(part_in[i], config_file)
    descr += arrow
    for i in xrange(part_out.shape[0]):
        descr += sb.pdg_to_name(part_out[i], config_file)
    if (clebsch_from_str(clebsch) != 1.0):
       descr = clebsch.replace(' ','') + '$\\times$(' + descr + ')'
    return descr

# Preambula configuration for plotting
execfile(os.path.dirname(os.path.abspath(__file__))+'/../../../python_scripts/common_plotting.py')

with open(input_file, 'r') as f:
    f.readline()
    smash_version = f.readline()
    f.readline()
    total_events = int(f.readline())
    f.readline()
    tend = float(f.readline())
    f.readline()
    tstart = float(f.readline())
    f.readline()
    reactions_string = f.readline().rstrip()
    first_reaction_descr = f.readline()
    bins_descr = first_reaction_descr.split(" bins")[0].split(', ')[-1]

contents = np.loadtxt(input_file, skiprows=10)

norm = total_events * (tend - tstart)

reactions_list = []
clebsch_factors_str = []
reactions_str = reactions_string.split('|')
for reaction_str in reactions_str:
    if ('x' in reaction_str):
        reaction_str_splitted = reaction_str.split('x')
        clebsch_factors_str.append(reaction_str_splitted[0])
        reaction_str = reaction_str_splitted[1]
    else:
        clebsch_factors_str.append('1')
    in_str, out_str = reaction_str.split(':')
    reac_in = np.sort(np.array([int(sb.name_to_pdg(x, config_file)) for x in in_str.split(',')]))
    reac_out = np.sort(np.array([int(sb.name_to_pdg(x, config_file)) for x in out_str.split(',')]))
    reactions_list.append([reac_in, reac_out])
react_num = len(reactions_list)
clebsch_factors_num = [clebsch_from_str(c) for c in clebsch_factors_str]

# colors and linestyle lists for plotting
col = cycle('brgkycm')
styles = cycle([':','-','--',])

# Set up subplot grid
gs = gridspec.GridSpec(40,4)
plt.subplot(gs[:11, :3])

descr_both = []
forward_backward_difference = []
plt.yscale('log')
largest_x = 0.0
smallest_x = 5.0

for i in xrange(react_num):
    binwidth = (contents[3*i][1] - contents[3*i][0])
    reaction_lhs = reactions_list[i][0]
    reaction_rhs = reactions_list[i][1]
    descr_both.append(get_descr(clebsch_factors_str[i], reaction_lhs, reaction_rhs, '$\\leftrightarrow$'))
    # Do not show points, where there are too little reactions in bin
    large_enough = (contents[3*i+1] > 5.0)
    M_inv = contents[3*i][large_enough]
    forward_rate = contents[3*i+1][large_enough]
    backward_rate = contents[3*i+2][large_enough]
    react_factor = clebsch_factors_num[i] / norm / binwidth
    color = next(col)
    style = next(styles)
    if (forward_rate > 0).any():
        plt.errorbar(M_inv, forward_rate * react_factor,
                       yerr=np.sqrt(forward_rate) * react_factor, color=color, fmt='>')
        if (M_inv[-1] > largest_x): largest_x = M_inv[-1]
        if (M_inv[0] < smallest_x): smallest_x = M_inv[0]
    else:
        print "WARN: no forward reactions for {} have been observed.".format(descr_both[i])
    if (backward_rate > 0).any():
        plt.errorbar(M_inv, backward_rate * react_factor,
                        yerr=np.sqrt(backward_rate) * react_factor, color=color, fmt='<', label = descr_both[i])
        if (M_inv[-1] > largest_x): largest_x = M_inv[-1]
        if (M_inv[0] < smallest_x): smallest_x = M_inv[0]
    else:
        print "WARN: no backward reactions for {} have been observed.".format(descr_both[i])

    # store difference between forward and backward reactions for second subplot (normalized to forward reactions)
    if (forward_rate > 0).any() and (backward_rate > 0).any():
        forward_backward_difference.append([M_inv,
        (forward_rate * react_factor - backward_rate * react_factor) / (forward_rate * react_factor),
        color, descr_both[i], style])

descr_isospinless = [isospin_form(react) for react in reactions_list]
descr_isospinless_unique = list(set(descr_isospinless))
in_which_isospin_group = [descr_isospinless_unique.index(descr_iso) for descr_iso in descr_isospinless]

n_isospin_groups = len(descr_isospinless_unique)
n_reactions_in_isospin_group = np.zeros(n_isospin_groups)
total_nscat_in_isospin_group = np.zeros(n_isospin_groups)
norm_isospin_group = np.zeros(react_num)
for i in xrange(react_num):
    i_group = in_which_isospin_group[i]
    n_reactions_in_isospin_group[i_group] += 2
    total_nscat_in_isospin_group[i_group] += contents[3*i+1].sum() * clebsch_factors_num[i]
    total_nscat_in_isospin_group[i_group] += contents[3*i+2].sum() * clebsch_factors_num[i]
for i in xrange(react_num):
    i_group = in_which_isospin_group[i]
    norm_isospin_group[i] = total_nscat_in_isospin_group[i_group]/n_reactions_in_isospin_group[i_group]


plt.minorticks_on()
smash_code_version, analysis_version = smash_version.split()
plt.figtext(0.8, 0.33, " %d events\n SMASH code:      %s\n SMASH analysis: %s\n"
            " time %.0f < t < %.0f fm/c" % \
            (total_events, smash_code_version, analysis_version, tstart, tend), \
            color = "gray", fontsize = 10)
plt.xticks(fontsize = 35)
plt.yticks(fontsize = 35)
plt.xlim(xmin = max(0.0, smallest_x - 0.2), xmax = largest_x + 0.2)
if (bins_descr == "M_inv"):
    plt.ylabel("$\\frac{dN_{react}}{dM_{inv} dt}$ [$GeV^{-1}$ $fm^{-1}$]", fontsize = 32)
elif (bins_descr == "|t|"):
    plt.ylabel("$\\frac{dN_{react}}{dt d\\tau}$ [$GeV^{-2}$ $fm^{-1}$]", fontsize = 32)
else:
    print "Unexpected bins description ", bins_descr


plt.subplot(gs[13:24, :3])
if (bins_descr == "M_inv"):
    plt.xlabel("$M_{inv}$, GeV", fontsize = 32)
    plt.ylabel("$\\frac{\\frac{dN_{react}^{\\qquad \\blacktriangleright}}{dM_{inv} dt} - \\frac{dN_{react}^{\\qquad \\blacktriangleleft}}{dM_{inv} dt}}{\\frac{dN_{react}^{\\qquad \\blacktriangleright}}{dM_{inv} dt}}$", fontsize = 32)
elif (bins_descr == "|t|"):
    plt.xlabel("$|t|$, [$GeV^{2}$]", fontsize = 32)
    plt.ylabel("$\\frac{\\frac{dN_{react}^{\\qquad \\blacktriangleright}}{dM_{inv} dt} - \\frac{dN_{react}^{\\qquad \\blacktriangleleft}}{dM_{inv} dt}}{\\frac{dN_{react}^{\\qquad \\blacktriangleright}}{dM_{inv} dt}}$", fontsize = 32)
else:
    print "Unexpected bins description ", bins_descr

for i in range(0,len(forward_backward_difference)):
    plt.plot(forward_backward_difference[i][0],
             forward_backward_difference[i][1],
             color = forward_backward_difference[i][2],
             label = forward_backward_difference[i][3],
             ls = forward_backward_difference[i][4],
             marker = 'o', lw = 3)

plt.axhline(y = 0, linestyle = '-', color = 'grey', zorder = 0, lw = 1)
plt.axhspan(-0.1, 0.1, color = 'grey', alpha = 0.2, zorder = 0, label = '10 %')
plt.xlim(xmin = max(0.0, smallest_x - 0.2), xmax = largest_x + 0.2)
plt.legend(prop={'size':25 if react_num < 10 else 19}, ncol = 1, frameon = 1, bbox_to_anchor=(1.365, 2.45))
plt.minorticks_on()
plt.xticks(fontsize = 35)
plt.yticks(fontsize = 35)

plt.subplot(gs[27:35, :])

def divide_ignore_zeros(a, b):
    """ Return c = a/b, but if some element b[i] = 0, then c[i] is assigned 0. """
    return np.divide(a, b, out = np.zeros_like(a), where = b!=0)

react_norm = divide_ignore_zeros(np.array(clebsch_factors_num), norm_isospin_group)
plt.errorbar(np.arange(react_num) - 0.05, contents[1::3].sum(axis=1)*react_norm, fmt='>',
             yerr=np.sqrt(contents[1::3].sum(axis=1))*react_norm, markersize=20.,
             capsize=20, elinewidth=4)
plt.errorbar(np.arange(react_num) + 0.05, contents[2::3].sum(axis=1)*react_norm, fmt='<',
             yerr=np.sqrt(contents[2::3].sum(axis=1))*react_norm, markersize=20.,
             capsize=20, elinewidth=4)

#print to file to store old_results
with open(results_file,"w") as f:
    f.write('channel' + '\t' + 'forward' + '\t' + 'backward' + '\n')
    for i in range(0,len(descr_both)):
        f.write(unicode(descr_both[i]) + '\t' + str((contents[1::3].sum(axis=1)*react_norm)[i]) + '\t' +
         str((contents[2::3].sum(axis=1)*react_norm)[i]) + '\n')
    f.write('# ' + smash_code_version)

ymax = max(max(contents[1::3].sum(axis=1)*react_norm),max(contents[2::3].sum(axis=1)*react_norm))

# plot reference data from previous SMASH version
if args.comp_prev_version:
    import comp_to_prev_version as cpv
    cpv.plot_previous_results('detailed_balance', setup, '/Nreact_by_Nisopingroup.txt', process_list = descr_both, ymax = ymax)

plt.plot(1,ymax * 4.0,linestyle = 'none', markersize = 20,
        color='black', marker = ">",label=str(smash_code_version))
plt.axhline(y = 1, linestyle = '-', color = 'grey', zorder = 0, lw = 1)
plt.xticks(range(react_num), descr_both, fontsize=15, rotation = 35)
plt.yticks([0,0.5,1,1.5,2,2.5,3,3.5,4], fontsize = 35)
plt.ylim(0.0, ymax * 1.5)
plt.xlim(-0.5, react_num-0.5)
plt.ylabel("$N_{react}/\langle N_{isospin\, group}\\rangle$", fontsize = 32)
plt.legend(prop={'size': 19}, loc='best', fontsize = 30, ncol = 2, frameon = 1)
plt.tight_layout(h_pad=100.0 )

plt.savefig(output_file)
