#!/usr/bin/python2

# angular.py
#
# This script reads SMASH binary output
# and prints out the angular distribution for processes
# which have the given PDG code(s) in initial state.
#
# Usage: ./angular.py [pdg code(s)] [output file(s)]
# Example: p+p: ./angular.py 2212 2212 data/0/collisions_binary.bin

import numpy as np
import sys       # system-specific functions
import os        # operating system interface
import argparse  # command line argument parser
import math      # math operations (sqrt, pi, ...)
from collections import defaultdict  # initialize dictionaries with certain type
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts')
import smash_basic_scripts as smash  # contains functions for reading SMASH binary

initial_pdgs = []                   # list of pdg codes (2 expected), read as input
theta_hist = defaultdict(int)       # histogram for distribution of theta (in various channels)
t_hist = defaultdict(int)           # histogram for distribution of t (in various channels)
theta_counts = defaultdict(int)
t_counts = defaultdict(int)
process_list = []                   # list of processes (= final-state channels)
thetabin = math.pi/30.              # bin size of theta histograms

parser = argparse.ArgumentParser()
# options and arguments
parser.add_argument("-v", "--verbose", action="store_true",
                    help="print additional output")
parser.add_argument("pdg1", help="PDG code of the first particle")
parser.add_argument("pdg2", help="PDG code of the second particle")
parser.add_argument("config_file", help="config file")
parser.add_argument("filename", nargs='+',
                    help="binary file(s) containing collision history")
args = parser.parse_args()

initial_pdgs.append(smash.pdg_to_name(int(args.pdg1), args.config_file))
initial_pdgs.append(smash.pdg_to_name(int(args.pdg2), args.config_file))
initial_pdgs.sort()

def abs3(p):
    "given a 4-vector p, determine the absolute value of its spatial 3-vector"
    return math.sqrt(p[1]**2 + p[2]**2 + p[3]**2)

def sqr4(p):
    "square a 4-vector p"
    return p[0]**2 - p[1]**2 - p[2]**2 -p[3]**2

def abs4(p):
    "determine the absolute value of a 4-vector p"
    return math.sqrt(sqr4(p))

def p_cm(srts,m1,m2):
    "get the center-of-mass momentum"
    s = srts**2
    m1_sqr = m1**2
    x = s + m1_sqr - m2**2
    return math.sqrt(x**2 * (0.25 / s) - m1_sqr)

def get_theta(p):
    "determine the polar angle theta of a 4-vector p"
    pabs = abs3(p)
    if (pabs>0):
      return math.acos(p[3]/pabs)
    else:
      return 0.

firstfile = True
num_elast = 0.  # number of elastic collisions
num_inel = 0.   # number of inelastic collsions
num_tot = 0.    # total number of collisions

# go through the list of input arguments
for name in args.filename:
    # if the argument is a file, read the contents
    if (os.path.isfile(name)):
        with smash.BinaryReader(name) as reader:
            # check the header for version information
            smashversion = reader.smash_version
            formatversion = reader.format_version
            if firstfile:
                smash_previous = smashversion
                format_previous = formatversion
                firstfile = False
            elif (smashversion != smash_previous or formatversion != format_previous):
                print "# Data from different versions detected!"
                print "# Above version from file", name
            # loop over all data blocks in the file
            for datablock in reader:
                # if interaction block, do the analysis
                if (datablock['type'] == 'i'):
                    # check that the initial state is what we want
                    initial = [smash.pdg_to_name(i) for i in datablock['incoming']['pdgid']]
                    if (sorted(initial) == initial_pdgs):
                        num_tot += 1.
                        final_pdgs = []
                        final_pdgs = [smash.pdg_to_name_generic(i) for i in datablock['outgoing']['pdgid']]
                        final_pdgs.sort()
                        process = "+".join(final_pdgs)
                        if process not in process_list:
                            process_list.append(process)
                        if (datablock['nout'] == 2):
                          # 2->2!
                          sigma_tot = datablock['total_cross_section']
                          # extract the 4-vectors of initial particles
                          p1 = datablock['incoming']['p'][0,:]
                          p2 = datablock['incoming']['p'][1,:]
                          p3 = datablock['outgoing']['p'][0,:]
                          p4 = datablock['outgoing']['p'][1,:]
                          # determine masses
                          m1 = abs4(p1)
                          m2 = abs4(p2)
                          m3 = abs4(p3)
                          m4 = abs4(p4)

                          # compute energy and momenta in center-of-mass
                          sqrts = abs4(p1 + p2)
                          p_i = p_cm(sqrts, m1, m2)
                          p_f = p_cm(sqrts, m3, m4)

                          # compute angle theta and mandelstam t
                          theta = get_theta(p3)
                          t_abs = -sqr4(p1-p3)

                          theta_index = int(theta / thetabin)
                          t_max = 4.*p_i*p_i   # maximum possible |t| in elastic collision
                          tbin = t_max / 50.
                          t_index = int(t_abs / tbin)
                          theta_hist[(process,theta_index)] += 1
                          theta_counts[theta_index] += 1
                          t_hist[(process,t_index)] += 1
                          t_counts[t_index] += 1

                          # clean up
                          del final_pdgs[:]

# sort the process list for consistency
process_list.sort()

### write out theta histogram
f = open('theta.dat', 'w+')
# print header
print >>f, "theta total",
for proc in process_list:
    print >> f, proc,
print >> f
# print data
for i in sorted(theta_counts):
  theta = (i+0.5)*thetabin    # use the center of the bin
  print >> f, theta, theta_counts[i]*sigma_tot/(num_tot*thetabin),
  for proc in process_list:
    print >> f, theta_hist[(proc,i)]*sigma_tot/(num_tot*thetabin),
  print >> f
print >>f, "#", smashversion

### write out t histogram
f = open('t.dat', 'w+')
# print header
print >>f, "t total",
for proc in process_list:
    print >> f, proc,
print >> f
# print data
for i in sorted(t_counts):
  t = (i+0.5)*tbin    # use the center of the bin
  print >> f , t, t_counts[i]*sigma_tot/(num_tot*tbin),
  for proc in process_list:
    print >>f , t_hist[(proc,i)]*sigma_tot/(num_tot*tbin),
  print >> f
print >>f, "#", smashversion
