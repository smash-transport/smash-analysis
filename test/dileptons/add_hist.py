import functools
import os
import numpy as np
from definitions import *  # IMPORT BINS, LINESTYLES AND CHANNELS

"""
Average histograms created by create_hist_from_event_output.py.
Plus, copy smash version number to parent directory.

"""


def numeric_compare(x, y):
    return int(x) - int(y)


def add(spectrum):
    n_folders = 0

    for i in sorted([f for f in os.listdir("data/") if not (f.startswith('.') or f == 'tabulations')], key=functools.cmp_to_key(numeric_compare)):
        # loop over folders, ignore hidden files
        if os.path.isdir("data/" + i):
            n_folders += 1

            filepath = "data/" + i + "/hist_raw_" + spectrum + ".txt"
            # print "Adding data folder %s ... " % filepath

            data = np.loadtxt(filepath)
            center = data[:,0]
            hist1 = data[:,1:]

            if n_folders == 1:
                hist_tot = hist1
            else:
                hist_tot += hist1

    hist_tot = hist_tot / n_folders

    np.savetxt("hist_raw_" + spectrum + ".txt", np.column_stack((center[:,None], hist_tot)))


def copy_version():
    # use first event
    data = np.genfromtxt("data/1/other.version.dat", dtype='str')
    with open("other.version.dat","w") as out:
        out.write(data[0] + " " + data[1] + "\n")


copy_version()

add("mass")
add("mass_rho")
add("mass_omega")
add("mass_0_800")
add("mass_800")

add("pt")
add("pt_rho")
add("pt_omega")
add("pt_0_150")
add("pt_150_470")
add("pt_470_700")
add("pt_700")

add("y")
add("y_rho")
add("y_omega")
add("y_0_150")
add("y_150_470")
add("y_470_700")
add("y_700")
