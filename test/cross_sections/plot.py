# -*- encoding: utf-8 -*-
import sys
import os
import argparse
import codecs
from collections import defaultdict, Counter
from itertools import cycle

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt


sys.path.append(
    os.path.dirname(
        os.path.abspath(__file__)) +
    '/../../python_scripts')
from common_plotting import smash_style, errorcontour, get_default_colors
from txt_io import load_table
from ordered_default_dict import OrderedDefaultDict
from collections import OrderedDict
from pdgs_from_config import strip_charge
import smash_basic_scripts as sb

reload(sys)
sys.setdefaultencoding('utf-8')

THIN_SPACE = u'\u2009'.encode('utf-8')

def parse_arguments():
    """Parse and return the command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Cross sections from SMASH")
    parser.add_argument("output_file", help="Output for plotting")
    parser.add_argument("--total_xs", help="Total cross sections from PDG", action="store")
    parser.add_argument("--elastic_xs", help="Elastic cross sections from PDG", action="store")
    parser.add_argument("--final_state", help="Only plot reactions with given ','-separated final states", type=unicode)
    parser.add_argument("--final_state_data", help="Plot given ','-separated experimental data files corresponding to the values given in `--final_state`.", type=str)
    parser.add_argument("system", help="Test setup for which plots are generated.", type = str)
    parser.add_argument("--comp_prev_version", required = False, action='store_true', help = "Plot comparison to previous SMASH version.")
    args = parser.parse_args()
    return args

def parse_state(particle_string):
    """Parse a string of particles delimited by '+' and return the count of each one."""
    return Counter(particle_string.split('+'))

def equal_state(left, right):
    """Determine whether the state `left` equals the state `right`.

    A state is a `Counter` of the particles in the state.
    'X' is intepreted as "any amount of other particles".
    """
    if left['X'] > 0 and right['X'] > 0:
        return True
    if right['X'] > 0:
        return all(left[p] >= right[p]
                   for p in right.keys() if p != 'X')
    if left['X'] > 0:
        return all(right[p] >= left[p]
                   for p in left.keys() if p != 'X')
    return (len(left) == len(right) and
            all(left[p] == right[p] for p in right.keys()))

def test_equal_state():
    left = parse_state('A+A+B+C')
    right = parse_state('A+B+C')
    assert not equal_state(left, right)
    assert not equal_state(right, left)
    right = parse_state('A+B+C+A')
    assert equal_state(left, right)
    assert equal_state(right, left)
    right = parse_state('A+C+X')
    assert equal_state(left, right)
    assert equal_state(right, left)
    right = parse_state('A+C+C+X')
    assert not equal_state(left, right)
    assert not equal_state(right, left)
    left = parse_state('B+D+X')
    assert equal_state(left, right)
    assert equal_state(right, left)

def dict_to_str(d):
    """Convert a dict to a human-readable string. Essential for debugging."""
    return '{' \
        + ', '.join('{}: {}'.format(k, v) for k, v in d.iteritems()) \
        + '}'

def list_to_str(l):
    """Convert a list to a human-readable string. Essential for debugging."""
    return '[' + ', '.join(str(i) for i in l) + ']'


if __name__ == '__main__':
    test_equal_state()

    args = parse_arguments()

    # (1) read reaction, initial masses, SMASH version and column names from header
    with codecs.open(args.input_file, 'r', encoding='utf-8') as f:
        reac = f.readline().rstrip().split()
        initial_masses = f.readline().rstrip().split()[1:]
        version = f.readline().rstrip().split()[1]
        colnames = f.readline().rstrip().split()[1:]

    fig = plt.figure()
    ax = plt.subplot(111)

    # (3) Plot experimental data
    if args.total_xs and os.path.isfile(args.total_xs):
        if os.path.isfile(args.total_xs):
            import comp_to_exp_data as ced
            ced.plot_cross_section_data(args.total_xs, initial_masses,
                                        total_xs = True)

        if args.elastic_xs and os.path.isfile(args.elastic_xs):
            if 'ced' not in sys.modules:
                import comp_to_exp_data as ced
            ced.plot_cross_section_data(args.elastic_xs, initial_masses,
                                        elastic_xs = True)
        # plot previous data here only for total

    # (4) Plot SMASH curves
    if args.final_state_data:
        final_state_data = args.final_state_data.split(',')
    else:
        final_state_data = None
    columns = np.loadtxt(args.input_file, skiprows=3, unpack=True)
    sqrts = columns[0]
    if args.final_state:
        final_state_names = args.final_state.split(',')
        final_states = [parse_state(s) for s in final_state_names]
        max_sigma = 0
        min_sqrts = np.amin(sqrts)
        max_sqrts = np.amax(sqrts)
    else:
        final_state_names = None
        final_states = None
        max_sigma = max(np.amax(columns[1:]), 80)
        min_sqrts = np.amin(sqrts)
        max_sqrts = np.amax(sqrts)
    if final_states and final_state_data:
        assert len(final_states) == len(final_state_data)

    nlabels = 1
    if len(columns) > 0:
        # Add up cross sections for cases where several final states are considered
        sigma_table = OrderedDefaultDict(lambda: np.zeros(len(sqrts)))
        plotted_data = defaultdict(bool)

        # Ignore the first column, which is sqrts
        # Skip every second column (contains error)
        for i in range(1, len(colnames), 2):
            if final_states:
                # Ignore the second and third column, which is sigma_tot and its error
                if i == 1:
                    continue
                # Sometimes we have "N*→Nπ", but we only care about the final state
                colname_split = colnames[i].split('→')
                assert len(colname_split) in (1, 2, 3), colname_split
                colname_without_res = colname_split[-1]
                state = parse_state(colname_without_res)
                for j, f in enumerate(final_states):
                    final_state_name = final_state_names[j]
                    if len(final_state_names) == 1:
                        final_state_name = 'exclusive'
                    final_state_name_err = final_state_name + '_err'
                    if equal_state(state, f):
                        if '→' not in colnames[i]:
                            # This is not a column with the initial resonance,
                            # so we have to be careful to account for "+X" final
                            # states.
                            sigma_table[final_state_name] += columns[i]
                            sigma_table[final_state_name_err] = np.sqrt(
                                sigma_table[final_state_name_err]**2 + columns[i + 1]**2)
                        else:
                            # This column has the initial resonance in its name.
                            # We still have to worry about "+X" final states to
                            # avoid blowing up the number of labels in the legend.
                            # We want to get rid of charges for the intermediate
                            # state.
                            lhs = colnames[i].split('→')[0]
                            lhs = '+'.join(map(strip_charge, lhs.split('+')))
                            colname = lhs + '→' + final_state_name
                            colname_err = colname + '_err'
                            sigma_table[colname] += columns[i]
                            sigma_table[colname_err] = np.sqrt(
                                sigma_table[colname_err]**2 + columns[i + 1]**2)
                        max_sigma = min(max(np.amax(sigma_table[final_state_name]),
                                            max_sigma), 80)
                        if final_state_data and not plotted_data[final_state_name]:
                            # plot experimental data
                            if 'ced' not in sys.modules:
                                import comp_to_exp_data as ced
                            max_sigma = ced.plot_cross_section_data(final_state_data[j],
                                                          initial_masses,
                                                          final_state_xs = True,
                                                          final_state_xs_name = final_state_name,
                                                          sqrts_range = [min_sqrts, max_sqrts],
                                                          max_sigma = max_sigma)

                            plotted_data[final_state_name] = True
            elif '→' not in colnames[i]:
                # If we plot everything, don't plot final states for each initial resonance,
                # the plot would get too crowded otherwise.
                sigma_table[colnames[i]] = columns[i]  # cross section
                sigma_table[colnames[i + 1]] = columns[i + 1]  # error
        total_sigma = sigma_table[u'total'].sum()
        if not total_sigma > 0.:
            # This only happens for plots of partial cross sections.
            # In such cases, we just assume the first cross section is the total.
            first = next(sigma_table.iterkeys())
            total_sigma = sigma_table[first].sum()
        assert total_sigma > 0.
        # Calculate contribution to total cross section
        sorted_sigma_table = {}
        for name, sigma in sigma_table.iteritems():
            if name.endswith('_err'):
                # We skip errors
                continue
            error = sigma_table[name + '_err']
            sorted_sigma_table[name] = (sigma, error, sigma.sum() / total_sigma)
        # Sort by contribution to total cross section
        sorted_sigma_table = OrderedDict(sorted(
             [(name, item)
              for name, item in sorted_sigma_table.iteritems()],
             key=lambda t: -t[1][2])
        )
        default_color = cycle(get_default_colors(plt.rcParams))
        plotted_lines = 0
        for name, t in sorted_sigma_table.iteritems():
            # Limit the number of plotted lines
            if plotted_lines > 100:
                break
            sigma, error, ratio = t
            # smoothing
            #from scipy.signal import savgol_filter
            #sigma = savgol_filter(sigma, 7, 3)
            # Don't show reactions in the legend if they contribute less than
            # 1% to the total cross section.
            if ratio > 0.01:
                if '→' in name:
                    label = name.split('→')[0]
                else:
                    label = name
                nlabels += 1
            else:
                label = '_nolegend_'
            label = label.replace('+', THIN_SPACE)
            color = next(default_color)
            errorcontour(sqrts, sigma, error, axis=ax, color=color,
                label=label)
            plotted_lines += 1
        if max_sigma == 0.:
            raise ValueError('nothing to plot')
        ax.set_ylim(0.0, max_sigma * 1.1)
        ax.set_xlim(min_sqrts, max_sqrts)

    # labels
    title = reac[1] + THIN_SPACE + reac[2]
    if final_state_names and len(final_state_names) == 1:
        title += '→' + final_state_names[0]
    title = title.replace('+', THIN_SPACE)
    ax.set_title(title)
    ax.set_xlabel(colnames[0])
    ax.set_ylabel("$\sigma$ [mb]")

    smash_style.set()

    # (4a) plot data from previous SMASH version
    # after applying the SMASH style! Otherwise this curve would be reformatted
    # plot only total_xs from prev. version if current one is plotted as well
    if args.comp_prev_version:
        import comp_to_prev_version as cpv
        cpv.plot_previous_results('cross_sections', str(args.system), '/sqrts_totalXSec.txt')

    # Shrink current axis
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.5, box.height])
    plt.legend(title=version, loc='center left', bbox_to_anchor=(1.04, 0.5),
               borderaxespad=0.)
    plt.figtext(0.8, 0.92, "SMASH analysis: %s" % \
                 (sb.analysis_version_string()), \
                 color = "gray", fontsize = 10)
    #if nlabels < 7:
    #    ncol = 1
    #elif nlabels < 13:
    #    ncol = 2
    #else:
    #    ncol = 3
    #ax.legend(loc="best", ncol=ncol)

    plt.savefig(args.output_file, bbox_inches='tight')
