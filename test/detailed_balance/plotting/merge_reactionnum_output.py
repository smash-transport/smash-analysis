#/usr/bin/python
# coding=UTF-8

import sys
import numpy as np
import argparse  # command line argument parser

parser = argparse.ArgumentParser()
parser.add_argument("output_file", type=str)
parser.add_argument("input_files", nargs='+',
                     help="text file(s) containing results of multiplicity vs time counting")
args = parser.parse_args()

def append_from_file(input_file, configuration, d):
    """
        Append dictionary d with with number of forward and reverse reactions,
        read from the file
    """
    with open(input_file, 'r') as f:
        _, smash_version, \
        _, total_events, \
        _, tend, \
        _, tstart, \
        _, reactions_str, \
        first_reaction_descr = \
        f.readline(), f.readline().rstrip(), \
        f.readline(), int(f.readline()), \
        f.readline(), float(f.readline()), \
        f.readline(), float(f.readline()), \
        f.readline(), f.readline().rstrip(), \
        f.readline()
    bins_descr = first_reaction_descr.split(" bins")[0].split(', ')[-1]
    contents = np.loadtxt(input_file, skiprows=10)
    to_test = (smash_version, total_events, tend, tstart)
    if (configuration == None):
        configuration = to_test
    elif (configuration != to_test):
        print "(Version, events, tend, tstart) are different in input files", \
              configuration, to_test
    reactions = reactions_str.split('|')
    for i in xrange(len(reactions)):
        r = reactions[i]
        # Get inverse reaction string
        if ('x' in r):
            m, a = r.split('x')
            m += 'x'
        else:
            m = ''
            a = r
        b, c = a.split(':')
        inverse_r = m + c + ':' + b
        # print r, inverse_r
        bins = contents[3*i]
        n_forward = contents[3*i + 1]
        n_reverse = contents[3*i + 2]
        if (r in d):
            d[r][1] += n_forward
            d[r][2] += n_reverse
        elif (inverse_r in d):
            d[inverse_r][2] += n_forward
            d[inverse_r][1] += n_reverse
        else:
            d[r] = [bins, n_forward, n_reverse]
    return configuration, bins_descr


conf = None
d = {}

for input_file in args.input_files:
    conf, bin_descr = append_from_file(input_file, conf, d)
# print conf, d

def reaction_order(s):
    isospinless = s.translate(None, '0+-⁺⁰⁻\^{}$')
    if 'x' in isospinless:
        isospinless = isospinless.split('x')[1]
    plus_count = s.count('⁺') + s.count('+')
    zero_count = s.count('⁰') + s.count('0')
    minus_count = s.count('⁻') + s.count('-')
    val = sum([ord(c) for c in isospinless])
    val += plus_count*3 + zero_count - minus_count
    return val

reactions_str = d.keys()
reactions_str.sort(key = lambda s: reaction_order(s))

reactions_string = '|'.join(reactions_str)
smash_version, events_per_file, tend, tstart = conf
total_events = events_per_file * len(args.input_files)

# Write output
with open(args.output_file, 'w') as f:
    f.write('# smash version\n')
    f.write('%s\n' % smash_version)
    f.write('# total number of events\n')
    f.write('%d\n' % total_events)
    f.write('# total time\n')
    f.write('%.2f\n' % tend)
    f.write('# time to start counting reactions\n')
    f.write('%.2f\n' % tstart)
    f.write('# list of all reaction analyzed\n')
    f.write('%s\n' % reactions_string)
with open(args.output_file, 'a') as f:
    for r in reactions_str:
        f.write('# reaction - %s, %s bins, forward and backward\n' % (r, bin_descr))
        np.savetxt(f, d[r][0], fmt = '%.3f', newline=' ')
        f.write('\n')
        np.savetxt(f, d[r][1], fmt = '%i', newline=' ')
        f.write('\n')
        np.savetxt(f, d[r][2], fmt = '%i', newline=' ')
        f.write('\n')
