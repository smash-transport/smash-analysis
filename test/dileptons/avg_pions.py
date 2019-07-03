import os
import numpy as np

import argparse

desc = """
        Returns pion event average.
        """

parser = argparse.ArgumentParser(description=desc)
parser.add_argument("particle_list", help="path to particle_lists.oscar")
args = parser.parse_args()

n_piz = 0  # number of pi^0
n_pip = 0  # number of pi^+
n_pim = 0  # number of pi^-

# print "Averaging pion for", args.particle_list

with open(args.particle_list) as f:

    for rawline in f:
        line = rawline.split()

        # particle info line
        if line[0][0] != "#":
            pdg = int(line[9])

            if pdg == 111:
                n_piz += 1
            elif pdg == 211:
                n_pip += 1
            elif pdg == -211:
                n_pim += 1

        # end-of-event marker
        if line[0] == "#" and line[1] == "event" and line[3] == "end":
            n_events = line[2]

# os.remove(args.particle_list)

with open(os.path.dirname(args.particle_list) + "/other.avg_pion.dat", "w") as out:
    # 1. line: average # of pi0
    out.write(str(float(n_piz) / float(n_events)) + "\n")
     # 2. line: average of (n_piz+n_pim)/2
    out.write(str(float(n_pip + n_pim) / (2.0 * float(n_events))) + "\n")
    # 3. line: # of events
    out.write(str(n_events) + "\n")
