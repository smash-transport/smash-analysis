import os
import numpy as np
import subprocess

"""
Average pion event average from different events.

"""


def numeric_compare(x, y):
    return int(x) - int(y)


def add_avg(particle):

    n_folders = 0

    for i in sorted([f for f in os.listdir("data/") if not f.startswith('.')], cmp=numeric_compare):
        # loop over folders, ignore hidden files
        if os.path.isdir("data/" + i):
            # print "Processing folder %s ... " % i
            n_folders += 1

            avgs1 = np.loadtxt("data/" + i + "/other.avg_" + particle + ".dat")

            if n_folders == 1:
                avgs = avgs1
            else:
                avgs = avgs + avgs1

    avgs_final = []
    avgs_final.append(avgs[0] / float(n_folders))
    avgs_final.append(avgs[1] / float(n_folders))
    avgs_final.append(avgs[2])  # n_events

    np.savetxt("other.avg_" + particle + ".dat", np.transpose(avgs_final))


add_avg("pion")
