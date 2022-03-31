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
    333: 17,                    # phi
    9922114: 18, 9922214: 18,  # N(2080)
    1218: 19,    2128: 19,     # N(2190)
    19922119: 20, 19922219: 20,  # N(2220)
    19932119: 21, 19932219: 21,  # N(2250)
    9962112: 22, 9962212: 22,  # N(1895)
    9922116: 23, 9922216: 23,  # N(2060)
    9972112: 24, 9972212: 24,  # N(2100)
    11212: 25,   12122:25,     # D(1900)
    "other": 26,
}

omega_channels = {
    21214: 0,    22124: 0,   # N(1700)
    42112: 1,    42212: 1,   # N(1710)
    31214: 2,    32124: 2,   # N(1720)
    9902114: 3,  9902214: 3,   # N(1875)
    9912114: 4,  9912214: 4,   # N(1900)
    9922114: 5,  9922214: 5,   # N(2080)
    1218: 6,     2128: 6,      # N(2190)
    9932114: 7 ,  9932214: 7,  # N(2120)
    9952112: 8,   9952212: 8,  # N(1880)
    9962112: 9,   9962212: 9,  # N(1895)
    9922116: 10,  9922216: 10, # N(2060)
    9972112: 11,  9972212: 11, # N(2100)
    "other": 12,
}

### CHANNEL LABELS ###

ch_list_main = [r'$\rho \rightarrow e^+e^-$',
                r'$\omega \rightarrow e^+e^-$',
                r'$\phi \rightarrow e^+e^-$',
                r'$\pi^0 \rightarrow \gamma e^+e^-$',
                r'$\eta \rightarrow \gamma e^+e^-$',
                r"$\eta' \rightarrow \gamma e^+e^-$",
                r'$\omega \rightarrow \pi^0 e^+e^-$',
                r'$\phi \rightarrow \pi^0 e^+e^-$',
                r'$\Delta^0\rightarrow n e^+e^-$',
                r'$\Delta^+ \rightarrow p e^+e^-$',
                r'$N^*(1520)\rightarrow N e^+e^-$',
                r'other']

ch_list_rho = [r'$\omega\rightarrow\rho\pi^0\rightarrow e^+e^-\pi^0$',
               r'$\pi^+\pi^-\rightarrow\rho\rightarrow e^+e^-$',
               r'$N^*(1520)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1535)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$\Delta^*(1620)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1650)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1675)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1680)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1700)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$\Delta^*(1700)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1710)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1720)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1875)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1900)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$\Delta^*(1905)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$\Delta^*(1950)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1990)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$\phi\rightarrow\rho\pi^0\rightarrow e^+e^-\pi^0$',
               r'$N^*(2080)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(2190)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(2220)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(2250)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(1895)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(2060)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$N^*(2100)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'$\Delta^*(1900)\rightarrow\rho N\rightarrow e^+e^-N$',
               r'other']

ch_list_omega = [r'$N^*(1700)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(1710)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(1720)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(1875)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(1900)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(2080)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(2190)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(2120)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(1880)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(1895)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(2060)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'$N^*(2100)\rightarrow\omega N\rightarrow e^+e^-N$',
                 r'other']


### LINESTYLES ###

line_style_main = ['b-', 'g-', 'r-', 'k--',
                   'c-', 'c--', 'g--', 'r--', 'y-', 'm--', 'b--' ,'k:']

# create (random) linestyle for origin plot
colors_o = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
linestyles_o = ['-', '--', '-.', ':']
nc = len(colors_o)
line_style_origin = []
for i in range(27):  # should be sufficent
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
