#!/usr/bin/python3

import os
from glob import glob
from txt_io import save_table, load_table
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict

"""
Collection of functions to easily compare latest SMASH results to those of
the previous version.

This script contains all the functions necessary to store, load and plot the
generated output files that are plotted in the different targets.
"""

def find_previous_version(observable, setup, filename):
    ''' Find existing SMASH data from previous runs, collect it for the given
        observable and sort it according to the version number in ascending
        order.
    '''

    labeled_ref_data = []
    pattern = os.path.dirname(os.path.abspath(__file__)) + '/../test/old_results/*/' + observable + '/' + setup + filename
    for path in glob(pattern):
        if not filename in path:
            continue
        label = [i for i in path.split('/') if 'SMASH' in i][0]

        ref_data = load_table(path)
        labeled_ref_data.append((label, ref_data))

    labeled_ref_data.sort(key=lambda x: x[0])
    return labeled_ref_data


def plot_previous_results(observable, setup, filename, color_list = [], pdg = 0,
                          process_list = [], plot_color = 'r', energy = 0.0,
                          plot_style = '-', ymax = 0.0, scaling_counter = 1.0):

    ''' This function is used to plot the data obtained with the preceding SMASH
        version in order to allow for an easy comparison. The only necessary
        arguments for the function call are observable, setup and filename. The
        corresponding data is loaded and plotted. For some observables
        additional arguments are necessary to properly compare to the latest
        results:
        • angular_distributions: color_list, process_list
        • detailed_balance: process_list, ymax
        • energy_scan (as function of sqrts): pdg, plot_color, plot_style
        • energy_scan (spectra): energy, plot_color, scaling_counter
    '''
    # (1) find reference data from previous version:
    labeled_reference_data = find_previous_version(observable, setup, filename)
    if labeled_reference_data == []:
        print('Warning: No reference data from previous version found for ' + str(observable) + ', ' + str(filename) + '. Comparison cannot be plotted.')
        return
    else:
        label = labeled_reference_data[-1][0]
        data = labeled_reference_data[-1][1]

        # (2) plot reference data
        if observable == 'angular_distributions':
            if 'theta' in filename:
                for i in range(0,len(data['theta'])):
                    s = np.sin(data['theta'][i])
                    for j in process_list:
                        if j in data:
                            data[j][i] = np.where(s > 0.0, data[j][i]/s, data[j][i])

                for i in process_list:
                    if i in data:
                        plt.plot(data['theta']*180/np.pi, data[i],
                                linestyle = '-', linewidth = 10, alpha = 0.15,
                                color = color_list.get(i),
                                zorder = 1)

            else:   # as function of t
                for i in process_list:
                    if i in data:
                        plt.plot(data['t'], data[i],
                                linestyle = '-', linewidth = 10, alpha = 0.15,
                                color = color_list.get(i),zorder = 1)

            #dummy for legend entry
            plt.plot(1,3,linestyle = '-', linewidth = 10,
                     alpha = 0.15, color='black',label=label)

        elif observable == 'cross_sections':
            plt.plot(data['sqrts'], data['x_sec'], linestyle = '-',
                     linewidth = 15, alpha = 0.2, color = "dimgrey",
                     zorder = 1, label = label + " (total)")

        elif observable == 'detailed_balance':
            channels = data['channel']
            previous_reactions_forward = data['forward'].astype(float)
            previous_reactions_backward = data['backward'].astype(float)

            for i,process in enumerate(channels):
                if process in process_list:
                    index = process_list.index(process)
                    plt.errorbar(index - 0.05, previous_reactions_forward[i],
                                fmt='4', markersize=22., capsize=20,
                                color = 'orange')
                    plt.errorbar(index + 0.05, previous_reactions_backward[i],
                                fmt='3', markersize=22., capsize=20,
                                color= 'deepskyblue')

            # dummy, just for legend entry
            plt.plot(1,ymax * 4.0,linestyle = 'none', markersize = 20,
                    color='black', marker = "4",label=label)

        elif observable == 'dileptons':
            plt.plot(data['m'], data['dN/dm'],
                    linestyle = '-', linewidth = 8, color = "dimgrey",
                    zorder = 1, alpha = 0.2, label = label + " (all)")

        elif observable == 'elastic_box':
            if 'geometric' in filename:
                plt.errorbar(data['x'], data['y'], yerr=data['y_error'],
                             fmt='s', ecolor = 'midnightblue', markeredgecolor= 'midnightblue', markersize = 15, capsize=10, label=label,
                             markerfacecolor= 'white', zorder = 1)
            elif 'stochastic' in filename:
                plt.errorbar(data['x'], data['y'], yerr=data['y_error'],
                             fmt='s', ecolor = 'maroon', markeredgecolor= 'maroon', markersize = 15, capsize=10, label=label,
                             markerfacecolor= 'white', zorder = 1)
            elif 'covariant' in filename:
                plt.errorbar(data['x'], data['y'], yerr=data['y_error'],
                             fmt='s', ecolor = 'forestgreen', markeredgecolor= 'forestgreen', markersize = 15, capsize=10, label=label,
                             markerfacecolor= 'white', zorder = 1)

        elif observable == 'energy_scan':
            if 'mtspectra' in filename:
                plt.plot(data['mt-m0'], data[str(energy)] * 10**scaling_counter,
                         color = plot_color, linestyle = '-', linewidth = 10,
                         zorder = 1, alpha = 0.15)
            elif 'yspectra' in filename:
                plt.plot(data['y'], data[str(energy)], color = plot_color,
                         linestyle = '-', linewidth = 10, zorder = 1, alpha = 0.15)
            elif 'ptspectra' in filename:
                plt.plot(data['pt'], data[str(energy)]* 10**scaling_counter, color = plot_color,
                         linestyle = '-', linewidth = 10, zorder = 1, alpha = 0.15)

            else:
                plt.plot(data['sqrt_s'], data[str(pdg)], plot_style,
                         color = plot_color, alpha = 0.15, linewidth = 10,
                         zorder = 1)
            return label

        elif observable == 'afterburner':
            if 'mtspectra' in filename:
                plt.plot(data['mt-m0'], data[str(energy)] * 10**scaling_counter,
                         color = plot_color, linestyle = '-', linewidth = 10,
                         zorder = 1, alpha = 0.15)
            elif 'yspectra' in filename:
                plt.plot(data['y'], data[str(energy)], color = plot_color,
                         linestyle = '-', linewidth = 10, zorder = 1, alpha = 0.15)
            elif 'ptspectra' in filename:
                plt.plot(data['pt'], data[str(energy)]* 10**scaling_counter, color = plot_color,
                         linestyle = '-', linewidth = 10, zorder = 1, alpha = 0.15)

            else:
                plt.plot(data['sqrt_s'], data[str(pdg)], plot_style,
                         color = plot_color, alpha = 0.15, linewidth = 10,
                         zorder = 1)
            return label


        elif observable == 'FOPI_pions':
            plt.plot(data['x'], data['y'], label=label, color = 'dimgrey',
                        linestyle = "-", alpha = 0.2, zorder=1, linewidth = 20)

        elif observable == 'densities':
            plt.plot(data['x'], data['y'], label=label, color = 'dimgrey',
                        linestyle = "-", alpha = 0.2, zorder=1, linewidth = 7)


        else:
            print('Error: Could not plot comparison to previous data. Unknown observable')
