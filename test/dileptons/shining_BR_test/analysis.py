import numpy as np
import os
import sys
import math
import argparse  # command line argument parser

desc = """
        Modified  and striped down dilepton analysis script from default analysis
        for testing whether the shining method reproduces the correct BR.
        Modified such that also pion invariant masses can be analysed.

        """

parser = argparse.ArgumentParser(description=desc)
parser.add_argument("data_file", help="file containing SMASH Oscar output")
parser.add_argument(
    "searched_pdg", help="positive pdg number of pair that is searched (works only for pions 211 or electons 11)", type=int)
args = parser.parse_args()

if args.searched_pdg != 11 and args.searched_pdg != 211:
    sys.exit("Unsupported PDG number entered.")

# only electron decays supported


class channel(object):
    rho = 0
    omega = 1
    phi = 2
    pi = 3
    eta = 4
    etaprime = 5
    omega_d = 6
    phi_d = 7
    delta0 = 8
    deltap = 9
    unknown = 10
    num = 11   # number of channels


def get_main_channel(n_out, pdg):
    if n_out == 2:  # 2-body decays
        if pdg == 113:
            return channel.rho
        elif pdg == 223:
            return channel.omega
        elif pdg == 333:
            return channel.phi
        else:
            return channel.unknown
    if n_out == 3:  # Dalitz decays
        if pdg == 111:
            return channel.pi
        elif pdg == 221:
            return channel.eta
        elif pdg == 331:
            return channel.etaprime
        elif pdg == 223:
            return channel.omega_d
        elif pdg == 333:
            return channel.phi_d
        elif pdg == 2214:
            return channel.deltap
        elif pdg == 2114:
            return channel.delta0
        else:
            return channel.unknown


def get_channel(dct, val):
    if val in dct:
        return dct[val]
    else:
        return len(dct.values())


# binning mass
bin_min = 0.0
bin_max = 2.0
nbins = 2000

bins_m = np.linspace(bin_min, bin_max, num=nbins + 1)
# you have to use 2 more histogram entries, because one is for counts below bin_min, the other for counts above bin_max
hist_mass = np.zeros((channel.num, nbins + 2))


n_folders = 0
num_events = 0
tmp_weight = 1.0  # standard weight (used for pions)

with open(args.data_file) as f:

    n_out = 0
    pdg = 0
    iterr = 0

    for rawline in f:
        line = rawline.split()

        #  header
        if line[0][0] == "#":

            if line[1] == "interaction":

                # check if iterator works corecctly
                if iterr != 0:
                    print "ERROR: Iterator fault"

                # grab weight, n_out and set iterator
                n_out = int(line[5])
                n_in = int(line[3])

                if args.searched_pdg == 11:
                    tmp_weight = float(line[9])
                iterr = n_out + n_in

            if line[1] == "event" and line[3] == "end":
                num_events = line[2]

        # particle line
        else:

                # grab particle info
            pdg = int(line[9])
            p = np.array([float(line[i]) for i in range(5, 9)])

            if iterr == n_out + 1 and n_in == 1 and not n_out == 1:  # decay, no wall crossing
                # determine channel
                in_part = get_main_channel(n_out, pdg)
            elif iterr < n_out + 1 and n_in == 1 and not n_out == 1:  # decay, no wall crossing
                # final state particles
                if pdg == args.searched_pdg:  # pi+ or e-
                    p_tmp = p
                elif pdg == -args.searched_pdg:  # pi- or e+
                    p_pi = [p_tmp, p]
                    # compute kinematic quantities
                    p_ges = p_pi[0] + p_pi[1]  # di-pion four-momentum
                    inv_mass = math.sqrt(
                        p_ges[0]**2 - p_ges[1]**2 - p_ges[2]**2 - p_ges[3]**2)
                    pt = math.sqrt(p_ges[1]**2 + p_ges[2]**2)
                    y = 0.5 * \
                        math.log((p_ges[0] + p_ges[3]) / (p_ges[0] - p_ges[3]))

                    # add to histograms
                    hist_mass[in_part, np.digitize(
                        [inv_mass], bins_m)] += tmp_weight

            iterr -= 1


def output(hist, bins, name):
    center = (bins[:-1] + bins[1:]) / 2.0
    # here the extra bin below bin_min and above bin_max is sliced out, so that the shapes of center and hist match again
    hist = hist[:, 1:-1] / float(num_events)
    with open("hist_" + name + ".txt", 'w') as result:
        np.savetxt(result, np.transpose(np.vstack([center, hist[:]])))


output(hist_mass, bins_m, "mass_" + str(args.searched_pdg))
