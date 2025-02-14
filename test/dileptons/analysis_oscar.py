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

# this is ordered as follows:
# - omega decay and pion production
# - nucleons and deltas up to 2 GeV in pole mass
# - light mesons up to 1.5 GeV in pole mass
# - remaining baryons
# within each, ordered by pole mass
rho_channels = {
       223: 0,                 # omega
       211: 1,      -211: 1,   # pion annihilation, no decay possible
      1214: 2,      2124: 2,   # N(1520)
     22112: 3,     22212: 3,   # N(1535)
      1212: 4,      2122: 4,   # D(1620)
     32112: 5,     32212: 5,   # N(1650)
      2116: 6,      2216: 6,   # N(1675)
     12116: 7,     12216: 7,   # N(1680)
     21214: 8,     22124: 8,   # N(1700)
     12114: 9,     12214: 9,   # D(1700)
     42112: 10,    42212: 10,  # N(1710)
     31214: 11,    32124: 11,  # N(1720)
   9902114: 12,  9902214: 12,  # N(1875)
   9902114: 13,  9902214: 13,  # N(1895)
   9912114: 14,  9912214: 13,  # N(1900)
      1216: 15,     2126: 14,  # D(1900)
      1216: 16,     2126: 15,  # D(1905)
      2118: 17,     2218: 16,  # D(1950)
   9902118: 18,  9902218: 18,  # N(1990)
       331: 19,                # η'
       333: 20,                # φ
     10223: 21,                # h₁(1170)
     10113: 22,    10213: 22,  # b₁(1235)
     20113: 23,    20213: 23,  # a₁(1260)
     10313: 24,    10323: 24,  # K₁(1270)
       225: 25,                # f₂
     20223: 26,                # f₁(1285)
    100111: 27,   100211: 26,  # π(1300)
       115: 28,      215: 27,  # a₂(1320)
     10221: 29,                # f₀(1370)
   9020221: 30,                # η(1405)
     20313: 31,    20323: 31,  # K₁(1400)
    100313: 32,   100323: 32,  # K*(1410)
    100223: 33,                # ω(1420)
       315: 34,      325: 34,  # K*₂(1430)
     10111: 35,    10211: 35,  # a₀(1450)
    100113: 36,   100213: 36,  # ρ(1450)
      3124: 37,                # Λ(1520)
     13124: 38,                # Λ(1690)
   9922116: 39,  9922216: 39,  # N(2060)
   9922114: 40,  9922214: 40,  # N(2080)
   9972112: 41,  9972212: 41,  # N(2100)
      1218: 42,     2128: 42,  # N(2190)
  19922119: 43, 19922219: 43,  # N(2220)
  19932119: 44, 19932219: 44,  # N(2250)
   "other": 45,
                             # higher N states ???
                             # other D states ???
  }

omega_channels = {
    21214: 0,    22124: 0,    # N(1700)
    42112: 1,    42212: 1,    # N(1710)
    31214: 2,    32124: 2,    # N(1720)
    9902114: 3,  9902214: 3,  # N(1875)
    9912114: 4,  9912214: 4,  # N(1900)
    9922114: 5,  9922214: 5,  # N(2080)
    1218: 6,     2128: 6,     # N(2190)
    9932114: 7 , 9932214: 7,  # N(2120)
    9952112: 8,  9952212: 8,  # N(1880)
    9962112: 9,  9962212: 9,  # N(1895)
    331: 10,                  # η'
    10113: 11,   10213: 11,   # b₁(1235)
    10313: 12,   10323: 12,   # K₁(1270)
    115: 13,     215:13,      # a₂(1320)
    20313: 14,   20323: 14,   # K₁(1400)
    100313: 15,  100323: 15,  # K*₂(1430)
    10111: 16,   10211: 16,   # a₀(1450)
    9922116: 17, 9922216: 17, # N(2060)
    9972112: 18, 9972212: 18, # N(2100)
    "other": 19,
}

phi_channels = {
  111: 0,   211:0,   -211:0,   # π ρ 
  113: 0,   213:0,   -213:0,   # π ρ 
  321: 1,   -321: 1,           # K⁺ K̅⁻ 
  311: 2,   -311: 2,           # K⁰ K̅⁰ 
  100333: 3,                   # φ(1680)
  9060225: 4,                  # f₂(2010)
  319: 5,   329: 5,            # K*₄(2045)
  9080225: 6,                  # f₂(2300)
  9090225: 7,                  # f₂(2340)
  9922114: 8,   9922214: 8,    # N(2080)
  9972112: 9,   9972212: 9,    # N(2100)
  9932114: 10,  9932214: 10,   # N(2120)
  1218:     11,    2128: 11,   # N(2190)
  19922119: 12, 19922219: 12,  # N(2220)
  19932119: 13, 19932219: 13,  # N(2250)
  "other": 14,
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
hist_mass_phi = np.zeros((max(phi_channels.values())+2, nbins+2))

# binning pt
bin_min = 0.0
bin_max = 2.0
nbins = 2000

bins_pt = np.linspace(bin_min, bin_max, num=nbins+1)
hist_pt = np.zeros((channel.num, nbins+2))
hist_pt_rho   = np.zeros((max(rho_channels.values())+2, nbins+2))
hist_pt_omega = np.zeros((max(omega_channels.values())+2, nbins+2))
hist_pt_phi = np.zeros((max(phi_channels.values())+2, nbins+2))

# binning rapidity
bin_min = -4.0
bin_max = 4.0
nbins = 2000

bins_rap = np.linspace(bin_min, bin_max, num=nbins+1)
hist_rap = np.zeros((channel.num, nbins+2))
hist_rap_rho   = np.zeros((max(rho_channels.values())+2, nbins+2))
hist_rap_omega = np.zeros((max(omega_channels.values())+2, nbins+2))
hist_rap_phi = np.zeros((max(phi_channels.values())+2, nbins+2))

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

print("lalalal",flush=True)
input()
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
                elif in_part == channel.phi:
                  # determine origin of phi meson
                  phi_ch = get_channel(phi_channels, parent)
                  hist_mass_phi[phi_ch, np.digitize([inv_mass], bins_m)] += tmp_weight
                  hist_pt_phi  [phi_ch, np.digitize([pt],      bins_pt)] += tmp_weight
                  hist_rap_phi [phi_ch, np.digitize([y],      bins_rap)] += tmp_weight

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
  output(hist_mass_phi, bins_m, "mass_phi")

  output(hist_pt,       bins_pt, "pt")
  output(hist_pt_rho,   bins_pt, "pt_rho")
  output(hist_pt_omega, bins_pt, "pt_omega")
  output(hist_pt_phi, bins_pt, "pt_phi")

  output(hist_rap,       bins_rap, "rapidity")
  output(hist_rap_rho,   bins_rap, "rapidity_rho")
  output(hist_rap_omega, bins_rap, "rapidity_omega")
  output(hist_rap_phi, bins_rap, "rapidity_phi")
