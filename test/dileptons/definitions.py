import numpy as np

"""
Definitions for dilepton analysis.

"""

### BINS ###

# binning mass
bin_min = 0.0
bin_max = 2.0
nbins_m = 150
bins_m = np.linspace(bin_min, bin_max, num=nbins_m + 1)

# binning pt
bin_min = 0.0
bin_max = 2.0
nbins_pt = 80
bins_pt = np.linspace(bin_min, bin_max, num=nbins_pt + 1)

# binning rapidity
bin_min = -4.0
bin_max = 4.0
nbins_y = 32
bins_y = np.linspace(bin_min, bin_max, num=nbins_y + 1)

# binning alpha
bin_min = -np.pi
bin_max = np.pi
nbins_alpha = 20  # very broad bins to save space
bins_alpha = np.linspace(bin_min, bin_max, num=nbins_alpha + 1)

centers_mass = (bins_m[:-1] + bins_m[1:]) / 2.0
centers_pt = (bins_pt[:-1] + bins_pt[1:]) / 2.0
centers_y = (bins_y[:-1] + bins_y[1:]) / 2.0
centers_alpha = (bins_alpha[:-1] + bins_alpha[1:]) / 2.0

hist_arr_dim = {"channel": 0,
                "mass": 1,
                "pt": 2,
                "y": 3,
                "alpha": 4}


### CHANNEL NUMBERING ###

n_channels = 12  # 11 known + 1 for potential unknown sources (should be zero)

direct_channels = {
    113:  0,  # rho
    223:  1,  # omega
    333:  2,  # phi
    "other": 11,  # unknown
}

dalitz_channels = {
    111:  3,  # pi
    221:  4,  # eta
    331:  5,  # etaprime
    223:  6,  # omega dalitz
    333:  7,  # phi dalitz
    2114:  8,  # delta 0
    2214:  9,  # delta p
    1214: 10,  2124: 10,  # N(1520)
    "other": 11,  # unknown
}

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

### CHANNEL LABELS ###

to_dil     = r'$\rightarrow e^+e^-$'
to_dil_N   = r'$\rightarrow N e^+e^-$'
to_dil_pi0 = r'$\rightarrow \pi^0 e^+e^-$'
to_dil_gamma = r'$\rightarrow \gamma e^+e^-$'

ch_list_main = [r'$\rho$',
                r'$\omega$' + to_dil,
                r'$\phi$' + to_dil,
                r'$\pi^0$' + to_dil_gamma,
                r'$\eta$' + to_dil_gamma,
                r"$\eta'$" + to_dil_gamma,
                r'$\omega$' + to_dil_pi0,
                r'$\phi$' + to_dil_pi0,
                r'$\Delta^0\rightarrow n e^+e^-$',
                r'$\Delta^+ \rightarrow p e^+e^-$',
                r'$N^*(1520)$' + to_dil_N,
                r'other']

to_rho_N = r'$\rightarrow N\rho$'
to_rho_pi0 = r'$\rightarrow \pi^0\rho$'
to_rho_eta = r'$\rightarrow \eta\rho$'
to_rho_omega = r'$\rightarrow \omega\rho$'
to_rho_K = r'$\rightarrow K\rho$'
to_rho_rho = r'$\rightarrow \rho\rho$'
ch_list_rho = [r'$\omega$' + to_rho_pi0,
               r'$\pi^+\pi^-\rightarrow\rho$',
               r'$N^*(1520)$' + to_rho_N, 
               r'$N^*(1535)$' + to_rho_N,
               r'$\Delta^*(1620)$' + to_rho_N,
               r'$N^*(1650)$' + to_rho_N,
               r'$N^*(1675)$' + to_rho_N,
               r'$N^*(1680)$' + to_rho_N,
               r'$N^*(1700)$' + to_rho_N,
               r'$\Delta^*(1700)$' + to_rho_N,
               r'$N^*(1710)$' + to_rho_N,
               r'$N^*(1720)$' + to_rho_N,
               r'$N^*(1875)$' + to_rho_N,
               r'$N^*(1895)$' + to_rho_N,
               r'$N^*(1900)$' + to_rho_N,
               r'$\Delta^*(1900)$' + to_rho_N,
               r'$\Delta^*(1905)$' + to_rho_N,
               r'$\Delta^*(1950)$' + to_rho_N,
               r'$N^*(1990)$' + to_rho_N,
               r'$\etaʹ\ \rightarrow\gamma\rho$',
               r'$\phi$' + to_rho_pi0,
               r'$h_1(1170)$' + to_rho_pi0,
               r'$b_1(1235)$' + to_rho_eta,
               r'$a_1(1260)$' + to_rho_pi0,
               r'$K_1(1270)$' + to_rho_K,
               r'$f_2$' + to_rho_rho,
               r'$f_1(1285)$' + to_rho_rho,
               r'$\pi(1300)$' + to_rho_pi0,
               r'$a_2(1320)\rightarrow(\pi^0,\omega)\rho$',
               r'$f_0(1370)$' + to_rho_rho,
               r'$\eta(1405)$' + to_rho_rho,
               r'$K_1(1400)$' + to_rho_K,
               r'$K^*(1410)\rightarrow K^*(892)\rho$',
               r'$\omega(1420)$' + to_rho_pi0,
               r'$K^*_2(1430)\rightarrow (K^*(892),K)\rho$',
               r'$a_0(1450)$' + to_rho_omega,
               r'$\rho(1450)\rightarrow(\eta,\rho)\rho$',
               r'$\Lambda(1520)\rightarrow\Sigma\rho$',
               r'$\Lambda(1690)\rightarrow\Sigma\rho$',
               r'$N^*(2060)$' + to_rho_N,
               r'$N^*(2080)$' + to_rho_N,
               r'$N^*(2100)$' + to_rho_N,
               r'$N^*(2190)$' + to_rho_N,
               r'$N^*(2220)$' + to_rho_N,
               r'$N^*(2250)$' + to_rho_N,
               r'other']

to_omega_N = r'$\rightarrow N\omega$'
to_omega_pi = r'$\rightarrow \pi\omega$'
to_omega_eta = r'$\rightarrow \eta\omega$'
to_omega_K = r'$\rightarrow K\omega$'
to_omega_gamma = r'$\rightarrow \gamma\omega$'
to_omega_rho = r'$\rightarrow \rho\omega$'
ch_list_omega = [r'$N^*(1700)$' + to_omega_N,
                 r'$N^*(1710)$' + to_omega_N,
                 r'$N^*(1720)$' + to_omega_N,
                 r'$N^*(1875)$' + to_omega_N,
                 r'$N^*(1900)$' + to_omega_N,
                 r'$N^*(2080)$' + to_omega_N,
                 r'$N^*(2190)$' + to_omega_N,
                 r'$N^*(2120)$' + to_omega_N,
                 r'$N^*(1880)$' + to_omega_N,
                 r'$N^*(1895)$' + to_omega_N,
                 r'$\etaʹ$' + to_omega_gamma,
                 r'$b_1(1235)$' + to_omega_pi,
                 r'$K_1(1270)$' + to_omega_K,
                 r'$a_2(1320)$' + to_omega_rho,
                 r'$K_1(1400)$' + to_omega_K,
                 r'$K^*_2(1430)$' + to_omega_K,
                 r'$a_0(1450)$' + to_omega_pi,
                 r'$N^*(2060)$' + to_omega_N,
                 r'$N^*(2100)$' + to_omega_N,
                 r'other']

to_phi = r'$\rightarrow \phi$'
to_phi_eta = r'$\rightarrow \eta\phi$'
to_phi_phi = r'$\rightarrow \phi\phi$'
to_phi_K892 = r'$\rightarrow K^*(892)\phi$'
to_phi_N = r'$\rightarrow N\phi$'
ch_list_phi = [r'$\pi\rho$' + to_phi, 
               r'$K^+K^-$' + to_phi, 
               r'$K^0\bar{K}^0$' + to_phi, 
               r'$\phi(1680)$' + to_phi_eta, 
               r'$f_2(2010)$' + to_phi_phi, 
               r'$K^*_4(2045)$' + to_phi_K892,    
               r'$f_2(2300)$' + to_phi_phi, 
               r'$f_2(2340)$' + to_phi_phi,             
               r'$N(2080)$' + to_phi_N,           
               r'$N(2100)$' + to_phi_N,           
               r'$N(2120)$' + to_phi_N,           
               r'$N(2190)$' + to_phi_N,           
               r'$N(2220)$' + to_phi_N,           
               r'$N(2250)$' + to_phi_N,
               r'other']

### LINESTYLES ###

line_style_main = ['b-', 'g-', 'r-', 'k--',
                   'c-', 'c--', 'g--', 'r--', 'y-', 'm--', 'b--' ,'k:']

# create (random) linestyle for origin plot
colors_o = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
linestyles_o = ['-', '--', '-.', ':']
nc = len(colors_o)
line_style_origin = []
for i in range(50):  # should be sufficent
    i = i % 26
    c = i % nc
    s = int(i / nc)
    line_style_origin.append(colors_o[c] + linestyles_o[s])


### PLOT STYLES ###

style_dict_mass = {'l_style': line_style_main,
                   'xlab': r'$m_{ee}\,[GeV]$',
                   'ylab': r'$dN/dm \,[1/GeV]$',
                   'x_min': 0.0,
                   'x_max': 1.2,
                   'y_min': 1E-8,
                   'y_max': 1E1}

style_dict_mass_origin = {'l_style': line_style_origin,
                          'xlab': r'$m_{ee}\,[GeV]$',
                          'ylab': r'$dN/dm \,[1/GeV]$',
                          'x_min': 0.0,
                          'x_max': 1.2,
                          'y_min': 1E-10,
                          'y_max': 1E-2}

style_dict_mass_w_data_pp_pNb = {'l_style': line_style_main,
                                 'xlab': r'$m_{ee}\,[GeV]$',
                                 'ylab': r'$d\sigma/dm\,[\mu b/GeV]$',
                                 'x_min': 0.0,
                                 'x_max': 1.2,
                                 'y_min': 1E-4,
                                 'y_max': 2E3}

style_dict_mass_w_data_CC_ArKCl = {'l_style': line_style_main,
                                   'xlab': r'$m_{ee}\,[GeV]$',
                                   'ylab': r'$1/N_{\pi^0}*dN/dm\,[1/GeV]$',
                                   'x_min': 0.0,
                                   'x_max': 1.2,
                                   'y_min': 1E-9,
                                   'y_max': 1E-2}

style_dict_pt = {'l_style': line_style_main,
                 'xlab': r'$p_T(e^+e^-)\,[GeV]$',
                 'ylab': r'$dN/dp_T \,[1/GeV]$',
                 'x_min': 0.0,
                 'x_max': 1.2,
                 'y_min': 1E-8,
                 'y_max': 1E-0}

style_dict_pt_origin = {'l_style': line_style_origin,
                        'xlab': r'$p_T(e^+e^-)\,[GeV]$',
                        'ylab': r'$dN/dp_T \,[1/GeV]$',
                        'x_min': 0.0,
                        'x_max': 1.2,
                        'y_min': 1E-10,
                        'y_max': 1E-2}

style_dict_pt_w_data = {'l_style': line_style_main,
                        'xlab': r'$p_T(e^+e^-)\,[GeV]$',
                        'ylab': r'$d\sigma/dp_T \,[\mu b/GeV]$',
                        'x_min': 0.0,
                        'x_max': 1.2,
                        'y_min': 1E-7,
                        'y_max': 1E2}

style_dict_y = {'l_style': line_style_main,
                'xlab': r'$y_{ee}$',
                'ylab': r'$dN/dy$',
                'x_min': -4.0,
                'x_max': 4.0,
                'y_min': 1E-8,
                'y_max': 1E0}

style_dict_y_origin = {'l_style': line_style_origin,
                       'xlab': r'$y_{ee}$',
                               'ylab': r'$dN/dy$',
                               'x_min': -4.0,
                               'x_max': 4.0,
                               'y_min': 1E-10,
                               'y_max': 1E-2}


style_dict_y_w_data = {'l_style': line_style_main,
                       'xlab': r'$y_{ee}$',
                       'ylab': r'$d\sigma/dy\,[\mu b]$',
                               'x_min': -4.0,
                               'x_max': 4.0,
                               'y_min': 1E-8,
                               'y_max': 1E1}
