import numpy as np
import os
import math
import argparse  # command line argument parser

parser = argparse.ArgumentParser()
parser.add_argument("data_file", help="file containing SMASH Oscar output")
parser.add_argument("acc_file", nargs='?', default="",
                    help="file containing acceptance matrix")
parser.add_argument("plep_min", type=float, nargs='?', default=0.,
                    help="minimum lepton momentum [GeV] (optional)")
parser.add_argument("plep_max", type=float, nargs='?', default=1000.,
                    help="maximum lepton momentum [GeV] (optional)")
parser.add_argument("res", type=int, nargs='?', default=3,
                    help="resolution mode (optional): 1=low, 2=medium, 3=high")
args = parser.parse_args()

# only electron decays supported
class channel(object):
    rho        = 0
    omega      = 1
    phi        = 2
    pi         = 3
    eta        = 4
    etaprime   = 5
    omega_d    = 6
    phi_d      = 7
    delta0     = 8
    deltap     = 9
    unknown    = 10
    num = 11   # number of channels

def get_main_channel(n_out, pdg):
  if n_out == 2:  # 2-body decays
    if   pdg == 113: return channel.rho
    elif pdg == 223: return channel.omega
    elif pdg == 333: return channel.phi
    else: return channel.unknown
  if n_out == 3:  # Dalitz decays
    if   pdg == 111:  return channel.pi
    elif pdg == 221:  return channel.eta
    elif pdg == 331:  return channel.etaprime
    elif pdg == 223:  return channel.omega_d
    elif pdg == 333:  return channel.phi_d
    elif pdg == 2214: return channel.deltap
    elif pdg == 2114: return channel.delta0
    else: return channel.unknown

rho_channels = {
      223: 0,                # omega
      211: 1,     -211: 1,   # pion annihilation, no decay possible
     1214: 2,     2124: 2,   # N(1520)
    22112: 3,    22212: 3,   # N(1535)
     1212: 4,     2122: 4,   # D(1620)
    32112: 5,    32212: 5,   # N(1650)
     2116: 6,     2216: 6,   # N(1675)
    12116: 7,    12216: 7,   # N(1680)
    21214: 8,    22124: 8,   # N(1700)
    12114: 9,    12214: 9,   # D(1700)
    42112: 10,   42212: 10,  # N(1710)
    31214: 11,   32124: 11,  # N(1720)
  9902114: 12, 9902214: 12,  # N(1875)
  9912114: 13, 9912214: 13,  # N(1900)
     1216: 14,    2126: 14,  # D(1905)
     2118: 15,    2218: 15,  # D(1950)
  9902118: 16, 9902218: 16,  # N(1990)
                             # higher N states ???
                             # other D states ???
  }

omega_channels = {
    21214: 0,    22124: 0,   # N(1700)
    42112: 1,    42212: 1,   # N(1710)
    31214: 2,    32124: 2,   # N(1720)
  9902114: 3,  9902214: 3,   # N(1875)
  9912114: 4,  9912214: 4,   # N(1900)
  9922114: 5,  9922214: 5,   # N(2080)
     1218: 6,     2128: 6,   # N(2190)
  }

def get_channel(dct, val):
  if val in dct:
    return dct[val]
  else:
    return max(dct.values())  # other

# binning mass
bin_min = 0.0
bin_max = 2.0
nbins = 2000

bins_m = np.linspace(bin_min, bin_max, num=nbins+1)
hist_mass = np.zeros((channel.num, nbins+2)) # you have to use 2 more histogram entries, because one is for counts below bin_min, the other for counts above bin_max
hist_mass_rho   = np.zeros((max(rho_channels.values())+2, nbins+2))  # maximum + 1, since 0 and +1, because of other
hist_mass_omega = np.zeros((max(omega_channels.values())+2, nbins+2))

# binning pt
bin_min = 0.0
bin_max = 2.0
nbins = 2000

bins_pt = np.linspace(bin_min, bin_max, num=nbins+1)
hist_pt = np.zeros((channel.num, nbins+2))
hist_pt_rho   = np.zeros((max(rho_channels.values())+2, nbins+2))
hist_pt_omega = np.zeros((max(omega_channels.values())+2, nbins+2))

# binning rapidity
bin_min = -4.0
bin_max = 4.0
nbins = 2000

bins_rap = np.linspace(bin_min, bin_max, num=nbins+1)
hist_rap = np.zeros((channel.num, nbins+2))
hist_rap_rho   = np.zeros((max(rho_channels.values())+2, nbins+2))
hist_rap_omega = np.zeros((max(omega_channels.values())+2, nbins+2))

n_folders = 0
num_events = 0


def get_abs_theta_phi(p):
  pabs = math.sqrt(p[1]**2 + p[2]**2 + p[3]**2)
  if (pabs > 0.):
    theta = math.acos(p[3]/pabs) * 180./math.pi
    phi   = math.atan2(p[2],p[1]) * 180./math.pi
    if (phi < 0.):
      phi += 360.
  else:
    theta = 0.
    phi   = 0.
  return pabs, theta, phi


# do HADES acceptance filtering
# input: list of four-vectors (electron and positron)
# note: the smeared four-vectors are being returned back to the caller
# output: acceptance probability
def HADES_filter(p):
  if not HADES_filter.enable:
    return 1.
  import HAFT  # Fortran-based HAFT code via f2py
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
  op_ang = math.acos( np.dot(p[0][1:],p[1][1:]) / (p1abs*p2abs) ) * 180./math.pi
  if (op_ang < 9.):
    return 0.
  # cut on single lepton momenta
  if (p1abs < args.plep_min or p1abs > args.plep_max or
      p2abs < args.plep_min or p2abs > args.plep_max):
    return 0.
  # determine acceptance probability
  acc = HAFT.haft_single.gethadesacceptance(3,p1abs,theta1,phi1,-2) * \
        HAFT.haft_single.gethadesacceptance(2,p2abs,theta2,phi2,-2)
  return acc
HADES_filter.enable = os.path.isfile(args.acc_file)  # enable acceptance filtering?
HADES_filter.init   = False                          # has the filter been initialized?


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
        if iterr != 0: print("ERROR: Iterator fault")

        # grab weight, n_out and set iterator
        n_out = int(line[5])
        tmp_weight = float(line[9])
        iterr = n_out + 1

      if line[1] == "event" and line[3] == "end":
        num_events = line[2]

    # particle line
    else:

        # grab particle info
        pdg = int(line[9])
        p = np.array([float(line[i]) for i in range(5,9)])

        if iterr == n_out + 1:
            # determine channel
            in_part = get_main_channel(n_out, pdg)
            # determine parent
            parent = int(line[17])
        elif iterr < n_out + 1:
            # final state particles
            if pdg == 11:  # electrons
              p_tmp = p
            elif pdg == -11:  # positrons
              p_lep = [p_tmp, p]
              acc = HADES_filter(p_lep)  # determine acceptance
              if acc > 0.:
                tmp_weight *= acc
                # compute kinematic quantities
                p_ges = p_lep[0] + p_lep[1]  # dilepton four-momentum
                inv_mass = math.sqrt(p_ges[0]**2 - p_ges[1]**2 - p_ges[2]**2 - p_ges[3]**2)
                pt = math.sqrt(p_ges[1]**2 + p_ges[2]**2)
                y = 0.5 * math.log((p_ges[0] + p_ges[3])/ (p_ges[0] - p_ges[3]))

                # add to histograms (mass, pt, rap)
                hist_mass[in_part, np.digitize([inv_mass], bins_m)] += tmp_weight
                hist_pt  [in_part, np.digitize([pt],      bins_pt)] += tmp_weight
                hist_rap [in_part, np.digitize([y],      bins_rap)] += tmp_weight

                if in_part == channel.rho:
                  # determine origin of rho meson
                  rho_ch = get_channel(rho_channels, parent)
                  hist_mass_rho[rho_ch, np.digitize([inv_mass], bins_m)] += tmp_weight
                  hist_pt_rho  [rho_ch, np.digitize([pt],      bins_pt)] += tmp_weight
                  hist_rap_rho [rho_ch, np.digitize([y],      bins_rap)] += tmp_weight
                elif in_part == channel.omega:
                  # determine origin of omega meson
                  omega_ch = get_channel(omega_channels, parent)
                  hist_mass_omega[omega_ch, np.digitize([inv_mass], bins_m)] += tmp_weight
                  hist_pt_omega  [omega_ch, np.digitize([pt],      bins_pt)] += tmp_weight
                  hist_rap_omega [omega_ch, np.digitize([y],      bins_rap)] += tmp_weight

        iterr -= 1


def output(hist, bins, name):
    center = (bins[:-1] + bins[1:]) / 2.0
    hist = hist[:,1:-1] / float(num_events) # here the extra bin below bin_min and above bin_max is sliced out, so that the shapes of center and hist match again
    with open("hist_" + name +".txt", 'w') as result:
      np.savetxt(result, np.transpose(np.vstack([center, hist[:]])))


if num_events==0:
  print("zero events found!")
else:
  #print num_events, "events"
  output(hist_mass,       bins_m, "mass")
  output(hist_mass_rho,   bins_m, "mass_rho")
  output(hist_mass_omega, bins_m, "mass_omega")

  output(hist_pt,       bins_pt, "pt")
  output(hist_pt_rho,   bins_pt, "pt_rho")
  output(hist_pt_omega, bins_pt, "pt_omega")

  output(hist_rap,       bins_rap, "rapidity")
  output(hist_rap_rho,   bins_rap, "rapidity_rho")
  output(hist_rap_omega, bins_rap, "rapidity_omega")
