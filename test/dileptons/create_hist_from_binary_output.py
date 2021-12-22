import numpy as np
import os
import sys
import math
import argparse  # command line argument parser
sys.path.append(os.path.dirname(
    os.path.abspath(__file__)) + '/../../python_scripts')
import smash_basic_scripts as sbs
from definitions import *  # IMPORT BINS, LINESTYLES AND CHANNELS

desc = """
        Analyse binary dilepton output and create histograms in inv. mass,
        pt, rapidity and azimutal angle, plus do HADES acceptance filtering.

        """

# NOTE: only di-eletrons supported
# NOTE: alpha = azimutal angle


### ARGUMENT PARSING ###
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("data_file", help="file containing SMASH dilepton output")
parser.add_argument("acc_file", nargs='?', default="",
                    help="file containing acceptance matrix")
parser.add_argument("plep_min", type=float, nargs='?', default=0.,
                    help="minimum lepton momentum [GeV] (optional)")
parser.add_argument("plep_max", type=float, nargs='?', default=1000.,
                    help="maximum lepton momentum [GeV] (optional)")
parser.add_argument("res", type=int, nargs='?', default=3,
                    help="resolution mode (optional): 1=low, 2=medium, 3=high")
args = parser.parse_args()


### DETERMINE CHANNEL (see definitions.py) ###
def get_channel(dct, val):
    if val in dct:
        return dct[val]
    else:
        return dct["other"]

def get_main_channel(n_out, pdg):
    if n_out == 2:  # 2-body decays
        return get_channel(direct_channels, pdg)
    if n_out == 3:  # Dalitz decays
        return get_channel(dalitz_channels, pdg)


### HADES ACCEPTANCE FILTERING ###
def get_abs_theta_phi(p):
    pabs = math.sqrt(p[1]**2 + p[2]**2 + p[3]**2)
    if (pabs > 0.):
        theta = math.acos(p[3] / pabs) * 180. / math.pi
        phi = math.atan2(p[2], p[1]) * 180. / math.pi
        if (phi < 0.):
            phi += 360.
    else:
        theta = 0.
        phi = 0.
    return pabs, theta, phi

# input: list of four-vectors (electron and positron)
# note: the smeared four-vectors are being returned back to the caller
# output: acceptance probability


def HADES_filter(p):
    if not HADES_filter.enable:
        return 1.
    import HAFT  # Fortran-based HAFT code via f2py3
    if not HADES_filter.init:
        HAFT.haft_single.setfilename(args.acc_file)
        HADES_filter.init = True
    # resolution smearing
    HAFT.haft_single.smearhadesmomentum(p[0], args.res, 3)   # 3=electron
    HAFT.haft_single.smearhadesmomentum(p[1], args.res, 2)   # 2=positron
    # get absolute momentum and angles
    p1abs, theta1, phi1 = get_abs_theta_phi(p[0])
    p2abs, theta2, phi2 = get_abs_theta_phi(p[1])
    # opening angle cut
    op_ang = math.acos(np.dot(p[0][1:], p[1][1:]) /
                       (p1abs * p2abs)) * 180. / math.pi
    if (op_ang < 9.):
        return 0.
    # cut on single lepton momenta
    if (p1abs < args.plep_min or p1abs > args.plep_max or
            p2abs < args.plep_min or p2abs > args.plep_max):
        return 0.
    # determine acceptance probability
    acc = HAFT.haft_single.gethadesacceptance(3, p1abs, theta1, phi1, -2) * \
        HAFT.haft_single.gethadesacceptance(2, p2abs, theta2, phi2, -2)
    return acc


HADES_filter.enable = os.path.isfile(
    args.acc_file)  # enable acceptance filtering?
# has the filter been initialized?
HADES_filter.init = False


### HISTOGRAM INIT ###

hist_mass       = np.zeros((nbins_m + 2, n_channels))
hist_mass_0_800 = np.zeros((nbins_m + 2, n_channels))
hist_mass_800   = np.zeros((nbins_m + 2, n_channels))
hist_mass_rho   = np.zeros((nbins_m + 2, max(rho_channels.values()) + 1))
hist_mass_omega = np.zeros((nbins_m + 2, max(omega_channels.values()) + 1))

hist_pt         = np.zeros((nbins_pt + 2, n_channels))
hist_pt_0_150   = np.zeros((nbins_pt + 2, n_channels))
hist_pt_150_470 = np.zeros((nbins_pt + 2, n_channels))
hist_pt_470_700 = np.zeros((nbins_pt + 2, n_channels))
hist_pt_700     = np.zeros((nbins_pt + 2, n_channels))
hist_pt_rho     = np.zeros((nbins_pt + 2, max(rho_channels.values()) + 1))
hist_pt_omega   = np.zeros((nbins_pt + 2, max(omega_channels.values()) + 1))

hist_y          = np.zeros((nbins_y  + 2, n_channels))
hist_y_0_150    = np.zeros((nbins_y  + 2, n_channels))
hist_y_150_470  = np.zeros((nbins_y  + 2, n_channels))
hist_y_470_700  = np.zeros((nbins_y  + 2, n_channels))
hist_y_700      = np.zeros((nbins_y  + 2, n_channels))
hist_y_rho      = np.zeros((nbins_y  + 2, max(rho_channels.values()) + 1))
hist_y_omega    = np.zeros((nbins_y  + 2, max(omega_channels.values()) + 1))


### ACTUAL ANALYSIS ###
unknown_ch = []
unknown_rho_ch = []
unknown_omega_ch = []


with sbs.BinaryReader(args.data_file) as reader:

    # save smash version number to extra data file
    smash_version = reader.smash_version
    with open("other.version.dat","w") as out:
        out.write("smash_version: ")
        out.write(smash_version + "\n")

    for block in reader:

        if block["type"] == 'f':
            num_events = block["nevent"]
        if block["type"] == 'p':
            sys.exit(
                "Error: Found particle block in dilepton binary collision output.")
        if block["type"] == 'i':
            if not (block["nin"] == 1 and block["nout"] != 1):  # not an decay
                print("Warning: Found a non-decay interaction block in Dilepton Output. (nin =", block["nin"], ", nout =", block["nout"], ")")
                continue

            shining_weight = block["total_cross_section"]
            decay_channel = get_main_channel(
                block["nout"], block["incoming"]["pdgid"][0])
            origin_pdg = block["incoming"]["PDG_mother1"][0]

            # Prevent double counting of omega dilepton Dalitz decay:
            # Currently we use a direct 3-body decay from the omega,
            # but we also have Dalitz like contributions coming from
            # omega decaying into rho pi, where the rho decays than
            # into dileptons. Those are blocked, for now.
            if decay_channel == 0 and origin_pdg == 223:
                continue

            if decay_channel == 10 and block["incoming"]["pdgid"][0] not in unknown_ch:
                unknown_ch.append(block["incoming"]["pdgid"][0])

            p_lep = [None] * 2
            # order is important for filtering
            for decay_product in block["outgoing"]:
                if decay_product["pdgid"] == 11:  # electron
                    p_lep[0] = decay_product["p"]
                elif decay_product["pdgid"] == -11:  # position
                    p_lep[1] = decay_product["p"]

            if p_lep[0] is None or p_lep[1] is None:
                raise ValueError("Found interaction without electron pair in decay products.")

            # make sure p_lep arrays satisfy fortran memory layout specifications
            p_lep = [np.require(p1, requirements=['F', 'A']) for p1 in p_lep]

            acc = HADES_filter(p_lep)  # determine acceptance
            if acc > 0.:
                shining_weight *= acc

                # compute kinematic quantities
                p_ges = p_lep[0] + p_lep[1]  # dilepton four-momentum
                inv_mass = math.sqrt(p_ges[0]**2 - p_ges[1]**2 - p_ges[2]**2 - p_ges[3]**2)
                pt = math.sqrt(p_ges[1]**2 + p_ges[2]**2)
                y = 0.5 * math.log((p_ges[0] + p_ges[3]) / (p_ges[0] - p_ges[3]))
                # alpha = math.atan2(p_ges[2], p_ges[1])  # azimutal angle (phi already used)
                p_ee = math.sqrt(p_ges[1]**2 + p_ges[2]**2 + p_ges[3]**2)

                # select bin and add shinning weight
                hist_mass[np.digitize([inv_mass], bins_m),  decay_channel] += shining_weight
                hist_pt  [np.digitize([pt],       bins_pt), decay_channel] += shining_weight
                hist_y   [np.digitize([y],        bins_y),  decay_channel] += shining_weight

                if 0.8 < p_ee: hist_mass_800   [np.digitize([inv_mass], bins_m), decay_channel] += shining_weight
                if p_ee < 0.8: hist_mass_0_800 [np.digitize([inv_mass], bins_m), decay_channel] += shining_weight

                if 0.00 < inv_mass < 0.15: hist_pt_0_150   [np.digitize([pt], bins_pt), decay_channel] += shining_weight
                if 0.15 < inv_mass < 0.47: hist_pt_150_470 [np.digitize([pt], bins_pt), decay_channel] += shining_weight
                if 0.47 < inv_mass < 0.70: hist_pt_470_700 [np.digitize([pt], bins_pt), decay_channel] += shining_weight
                if 0.70 < inv_mass:        hist_pt_700     [np.digitize([pt], bins_pt), decay_channel] += shining_weight

                if 0.00 < inv_mass < 0.15: hist_y_0_150    [np.digitize([y],  bins_y),  decay_channel] += shining_weight
                if 0.15 < inv_mass < 0.47: hist_y_150_470  [np.digitize([y],  bins_y),  decay_channel] += shining_weight
                if 0.47 < inv_mass < 0.70: hist_y_470_700  [np.digitize([y],  bins_y),  decay_channel] += shining_weight
                if 0.70 < inv_mass:        hist_y_700      [np.digitize([y],  bins_y),  decay_channel] += shining_weight

                # determine origin of rho and omega meson
                if decay_channel == 0:  # rho
                    rho_ch = get_channel(rho_channels, origin_pdg)
                    if rho_ch == 0:
                        print("Warning: Potential double counting. Found rho dilepton decay originating from omega.")
                    elif rho_ch == rho_channels["other"] and origin_pdg not in unknown_rho_ch:
                        unknown_rho_ch.append(origin_pdg)
                    hist_mass_rho [np.digitize([inv_mass], bins_m),  rho_ch] += shining_weight
                    hist_pt_rho   [np.digitize([pt],       bins_pt), rho_ch] += shining_weight
                    hist_y_rho    [np.digitize([y],        bins_y),  rho_ch] += shining_weight
                elif decay_channel == 1:  # omega
                    omega_ch = get_channel(omega_channels, origin_pdg)
                    if omega_ch == omega_channels["other"] and origin_pdg not in unknown_omega_ch:
                        unknown_omega_ch.append(origin_pdg)
                    hist_mass_omega [np.digitize([inv_mass], bins_m),  omega_ch] += shining_weight
                    hist_pt_omega   [np.digitize([pt],       bins_pt), omega_ch] += shining_weight
                    hist_y_omega    [np.digitize([y],        bins_y),  omega_ch] += shining_weight

if unknown_ch != []:       print("Warning: Unknown dilepton decay(s) found! -->", unknown_ch)
if unknown_rho_ch != []:   print("Warning: Unknown dilepton decay origin(s) for rho! -->", unknown_rho_ch)
if unknown_omega_ch != []: print("Warning: Unknown dilepton decay origin(s) for omega! -->", unknown_omega_ch)

num_events = int(num_events) + 1  # FIXME: event counting starts at zero in binary output

### OUTPUT ###
def output(hist, bins, name):
    center = (bins[:-1] + bins[1:]) / 2.0
    hist = hist[1:-1,:] / float(num_events) # here the extra bin below bin_min and above bin_max is sliced out, so that the shapes of center and hist match again
    np.savetxt("hist_" + name +".txt", np.column_stack((center[:,None], hist)))


output(hist_mass,       bins_m, "mass")
output(hist_mass_rho,   bins_m, "mass_rho")
output(hist_mass_omega, bins_m, "mass_omega")
output(hist_mass_0_800, bins_m, "mass_0_800")
output(hist_mass_800,   bins_m, "mass_800")

output(hist_pt,         bins_pt, "pt")
output(hist_pt_rho,     bins_pt, "pt_rho")
output(hist_pt_omega,   bins_pt, "pt_omega")
output(hist_pt_0_150,   bins_pt, "pt_0_150")
output(hist_pt_150_470, bins_pt, "pt_150_470")
output(hist_pt_470_700, bins_pt, "pt_470_700")
output(hist_pt_700,     bins_pt, "pt_700")

output(hist_y,         bins_y, "y")
output(hist_y_rho,     bins_y, "y_rho")
output(hist_y_omega,   bins_y, "y_omega")
output(hist_y_0_150,   bins_y, "y_0_150")
output(hist_y_150_470, bins_y, "y_150_470")
output(hist_y_470_700, bins_y, "y_470_700")
output(hist_y_700,     bins_y, "y_700")


# delete Dilepton Output
os.remove(args.data_file)
