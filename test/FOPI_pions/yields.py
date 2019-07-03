#!/usr/bin/python2
#
# This script determines pion yields from the binary particle output of SMASH.
#
# Usage: python yields.py > pion_yields.txt

import sys
import os
import os.path as path
from math import sqrt
from collections import defaultdict  # initialize dictionaries with certain type
import yaml
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts')
import smash_basic_scripts as smash  # contains functions for reading SMASH binary

# All particles that will be counted and written to the output file.
# If you change this, you have to change the plotting scripts as well!
pdglist_print = [211, 111, -211]

firstfile = True

class Average:
    def __init__(self):
        self.avg = 0.
        self.n = 0
        self.v = 0.

    def add(self, x):
        self.n += 1
        prev_avg = self.avg
        self.avg += (x - prev_avg) / self.n
        self.v += (x - prev_avg) * (x - self.avg)

    def var(self):
        if self.n < 2:
            return 0.
        return self.v / (self.n - 1)

    def err(self):
        if self.n == 0:
            return 0.
        return sqrt(self.var() / self.n)

    def __repr__(self):
        return self.avg

# make sure the implementation of Average is correct
a = Average()
for i in range(1, 6):
    a.add(i)
assert(a.avg == 3.0)
assert(a.var() == 2.5)

configpath = os.path.join(
    os.path.dirname('.'), '0.4', 'data', '1', 'config.yaml')
with open(configpath, 'r') as f:
    d = yaml.load(f)
ntestparticles = d['General']['Testparticles']

for edir in sorted(os.listdir(".")):    # loop over energies
    if path.isdir(edir) and (edir != "CMakeFiles"):
        datadir = path.join(path.abspath(edir), "data")
        count_with_error = defaultdict(Average)  # total count of processes with the initial pdgs
        pt_total = defaultdict(Average)  # average transverse momentum
        nevent = 0
        for rdir in os.listdir(path.abspath(datadir)):  # loop over runs
            rdir = path.join(path.abspath(datadir), rdir)
            if path.isdir(rdir):
                #print rdir
                name = path.join(rdir, "particles_binary.bin")
                # if the argument is a file, read the contents
                if (path.isfile(name)):
                    with smash.BinaryReader(name) as reader:
                        # check the header for version information
                        smashversion = reader.smash_version
                        formatversion = reader.format_version
                        if firstfile:
                            print "#version", smashversion, "format", formatversion, "analysis", smash.analysis_version_string()
                            firstfile = False
                        # loop over all data blocks in the file
                        for datablock in reader:
                            if datablock['type'] == 'p':
                                nevent += 1
                                # count particles and calculate transverse momentum
                                pt_event = defaultdict(Average)
                                count = defaultdict(int)
                                for particle in datablock['part']:
                                    i = particle['pdgid']
                                    if i not in pdglist_print:
                                        continue
                                    count[i] += 1
                                    px = particle['p'][1]
                                    py = particle['p'][2]
                                    pt_event[i].add(sqrt(px**2 + py**2))
                                for i in pdglist_print:
                                    count_with_error[i].add(count[i])
                                # add average pt of event to average of average pt of all events
                                for i in pdglist_print:
                                    # If no particle has been observed, it should not affect the
                                    # average pt.
                                    if pt_event[i].n > 0:
                                        pt_total[i].add(pt_event[i].avg)

        #print >> sys.stderr, edir, nevent       # debug output
        if nevent > 0:
            # print out energy, Npi+, Npi0, Npi-, Nevent, Ntestparticles, ptpi+, ptpi0, ptpi-, erorrs
            print edir,
            for i in pdglist_print:
                print count_with_error[i].avg, count_with_error[i].err(),
            print nevent, ntestparticles,
            for i in pdglist_print:
                print pt_total[i].avg, pt_total[i].err(),
            print
