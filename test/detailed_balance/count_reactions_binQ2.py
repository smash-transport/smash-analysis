import sys
import os  # operating system interface
import numpy as np
import argparse  # command line argument parser
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts')
import smash_basic_scripts as sb
from collections import defaultdict  # initialize dictionaries with certain type
import yaml


class Hist:
    """ Implements simple histogram for binning reactions. """

    def __init__(self, xmin, xmax, nbins):
        self.x_min = xmin
        self.x_max = xmax
        self.nbins = nbins
        self.dx_inv = float(self.nbins)/(self.x_max - self.x_min)
        self.hist = [0]*self.nbins

    def bin_centers(self):
        w = (self.x_max - self.x_min)/self.nbins
        return np.linspace(self.x_min + 0.5*w, self.x_max - 0.5*w, self.nbins)

    def fill(self, x):
        i = np.trunc((x - self.x_min) * self.dx_inv).astype(int)
        if (i >= 0 and i < self.nbins):
            self.hist[i] += 1
        else:
            print x, " is out of histogram range..."

def reaction_from_string(s, config_file):
    """
    Given a string s of the format n x a,b,...:c,d,... returns tuple ((a,b,),(c,d,)).
    Here : separates left and right hand sides, number before x and x itself are ignored,
    they are coefficients to compensate Clebsch-Gordan and symmetry factors.

    config file is only needed to get pdgcodes from particle names.
    """

    in_str, out_str = s.split(':')
    if 'x' in in_str:
        in_str = in_str.split('x')[1]
    reac_in  = [int(sb.name_to_pdg(x, config_file)) for x in  in_str.split(',')]
    reac_out = [int(sb.name_to_pdg(x, config_file)) for x in out_str.split(',')]
    reac_in.sort()
    reac_out.sort()
    return ((tuple(reac_in), tuple(reac_out)))

def tuple_to_name(tup, config_file):
    """ Reverse function to reaction_from_string, returns name given a tuple. """

    part_in, part_out = tup
    s  = ','.join([sb.pdg_to_name(x, config_file) for x in part_in]) + ':'
    s += ','.join([sb.pdg_to_name(x, config_file) for x in part_out])
    return s

def inv_reaction(r):
    """ Exchanges in and out in the tuple representing reaction. """
    return (r[1], r[0])


def count_reactions(files_to_analyze, tstart, only_count_reactions = None, binning_in = 'Q2'):
    """
    Counts, how many reactions of every possible type (binned by invariant mass)
    happen in files_to_analyze starting from time tstart.

    One may limit the list of considered reactions to only_count_reactions.
    """

    intcounter = 0  # any interactions counter
    event_num = 0   # event counter

    assert((binning_in == 'Q2') or (binning_in == 't'))
    def right_hist():
        if (binning_in == 'Q2'):
            return Hist(0., 5., 100)
        if (binning_in == 't'):
            return Hist(0., 30., 60)
    reactions = defaultdict(right_hist)

    for file_to_analyze in files_to_analyze:
        # Read input file, add Q^2 of given reactions to list: forward and backward
        bfile = open(file_to_analyze, 'rb')
        with sb.BinaryReader(file_to_analyze) as reader:
            smash_version = reader.smash_version
            #format_version = reader.format_version

            for block in reader:
                if (block['type'] == 'f'):  # end of event
                    event_num += 1
                if (block['type'] == 'i'):  # interaction
                    intcounter += 1
                    #  if (intcounter % 10000 == 0): print "interaction ", intcounter
                    time = sb.get_block_time(block)
                    if (time < tstart): continue
                    block_pdgin  = np.sort(block['incoming']['pdgid'])
                    block_pdgout = np.sort(block['outgoing']['pdgid'])
                    # ignore wall crossings
                    if ((block_pdgin.size == 1) and (block_pdgout.size == 1)): continue
                    # computing Mandelstam t requires 2->2 reaction
                    if ((binning_in == 't') and
                       (not ((block_pdgin.size == 2) and (block_pdgout.size == 2)))): continue
                    reaction = (tuple(block_pdgin), tuple(block_pdgout))
                    if (only_count_reactions and (not reaction in only_count_reactions)):
                        continue
                    if (binning_in == 'Q2'):
                        bin_variable = np.sqrt(sb.reaction_Q2(block))
                    if (binning_in == 't'):
                        bin_variable =  np.abs(sb.reaction_mandelstam_t(block))
                    reactions[reaction].fill(bin_variable)
    return (reactions, event_num, smash_version)

def find_n_most_detbal_violating(reactions, n = 6):
    deviations = dict()
    for r in reactions:
        n_forward = sum(reactions[r].hist)
        inverse_r = inv_reaction(r)
        if inverse_r in reactions:
            n_reverse = sum(reactions[inverse_r].hist)
        else:
            n_reverse = 0
        if ((r not in deviations) and (inverse_r not in deviations)):
            deviations[r] = np.abs(np.sqrt(n_forward) - np.sqrt(n_reverse))
    reactions_by_deviation = deviations.keys()
    reactions_by_deviation.sort(key = lambda(x): deviations[x])
    n = min(n, len(reactions_by_deviation))
    most_deviating_reactions = reactions_by_deviation[-n:]
    for r in most_deviating_reactions:
        print r, sum(reactions[r].hist), sum(reactions[(r[1], r[0])].hist), deviations[r]
    return most_deviating_reactions

def count_reactions_and_printout(binvar = 'Q2'):
    assert(binvar == 'Q2' or binvar == 't')
    parser = argparse.ArgumentParser()
    # options and arguments
    parser.add_argument("output_file", type=str)
    parser.add_argument("reactions_string", type=str)
    parser.add_argument("tstart", type=float)
    parser.add_argument("config_file", help="config file")
    parser.add_argument("files_to_analyze", nargs='+',
                        help="binary file(s) containing collision history")
    args = parser.parse_args()
    with open(args.config_file, 'r') as f: d = yaml.load(f)

    if (args.reactions_string == "n_most_violating"):
        reactions, event_num, smash_version = count_reactions(args.files_to_analyze, args.tstart, binning_in = binvar)
        most_deviating_reactions = find_n_most_detbal_violating(reactions)
        print most_deviating_reactions
        reaction_names = [tuple_to_name(x, args.config_file) for x in most_deviating_reactions]
    else:
        reaction_names = args.reactions_string.split('|')
        only_count_reactions = []
        for reaction_name in reaction_names:
            r = reaction_from_string(reaction_name, args.config_file)
            only_count_reactions.append(r)
            only_count_reactions.append(inv_reaction(r))
        reactions, event_num, smash_version = count_reactions(args.files_to_analyze, args.tstart,
                                               only_count_reactions = only_count_reactions, binning_in = binvar)

    # Write output
    with open(args.output_file, 'w') as f:
        f.write('# smash and analysis version\n')
        f.write('%s %s\n' % (smash_version, sb.analysis_version_string()))
        f.write('# total number events\n')
        f.write('%d\n' % event_num)
        f.write('# total time\n')
        f.write('%.2f\n' % d['General']['End_Time'])
        f.write('# time to start counting reactions\n')
        f.write('%.2f\n' % args.tstart)
        f.write('# list of all reaction analyzed\n')
        f.write('%s\n' % '|'.join(reaction_names))
    with open(args.output_file, 'a') as f:
        for rname in reaction_names:
            r = reaction_from_string(rname, args.config_file)
            binvar_name = {'Q2' : 'M_inv', 't' : '|t|'}
            f.write('# reaction - %s, %s bins, forward and backward\n' % (rname, binvar_name[binvar]))
            x = reactions[r].bin_centers()
            yf = reactions[r].hist
            yb = reactions[inv_reaction(r)].hist
            np.savetxt(f, x, fmt = '%.3f', newline=' ')
            f.write('\n')
            np.savetxt(f, yf, fmt = '%i', newline=' ')
            f.write('\n')
            np.savetxt(f, yb, fmt = '%i', newline=' ')
            f.write('\n')

if __name__ == '__main__':
    count_reactions_and_printout()
