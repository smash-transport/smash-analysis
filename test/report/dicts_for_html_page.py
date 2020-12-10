# coding=UTF-8

'''This script contains all known deviations of SMASH results. The dictionary is
sorted by the target first and the collision system second. It can be called
from the generate_html.py script. '''

target_dict = {
    'angular_distributions': 'Angular Distributions',
    'cross_sections': 'Cross Sections',
    'detailed_balance': 'Detailed Balance',
    'dileptons': 'Dileptons',
    'elastic_box': 'Scattering Rates',
    'FOPI_pions': 'FOPI Pions',
    'pp_collisions': 'PP Collisions',
    'spectra': 'Spectra',
    'densities': 'Densities',
    'energy_scan': 'Energy Scan',
    'total_multiplicity': 'Total Multiplicity',
    'midrapidity_yield' : 'Midrapidity Yield',
    'meanmt' : 'Mean m<sub>T</sub>',
    'meanpt' : 'Mean p<sub>T</sub>',
    'mtspectra' : 'm<sub>T</sub> Spectra',
    'yspectra' : 'Rapidity Spectra',
}

comments = {'cross_sections': {
        'kplus_neutron': 'Uses GiBUU parametrization. Does not describe exclusive cross section well, but it is very small and not important. See strangeness paper.',
        'kplus_proton': 'Exclusive Deltas are overestimated at high energies, probably because we don\'t produce enough pions.',
        'kminus_neutron': '-',
        'kminus_proton': 'Lambda(1520) peak could be higher, but it is tightly constrained by PDG branching ratios. At high energies, there are no resonances.',
        'neutron_proton': '-',
        'proton_proton': '-',
        'antiproton_proton': 'Failed strings correspond to processes in which pythia does not find a proper final state for an underlying string process. In this case SMASH falls back to an elastic scattering of the incoming particles.',
        'piplus_proton': 'Hyperon production discussed in <a href="https://arxiv.org/pdf/1809.03828.pdf"> strangeness paper</a>.',
        'piplus_piminus': '-',
        'piminus_proton': 'Hyperon production discussed in <a href="https://arxiv.org/pdf/1809.03828.pdf"> strangeness paper</a>.',
    }, 'elastic_box': {
        'scatrate_vs_N': '-',
        'scatrate_vs_Ntest': '-',
        'scatrate_vs_T': '-',
        'scatrate_vs_V': '-',
        'scatrate_vs_dt': '-',
        'scatrate_vs_sigma': '-',
    }, 'FOPI_pions': 'See <a href="http://inspirehep.net/record/1471536"> SMASH paper</a> (Fig. 25) for further discussion.',
       'detailed_balance': {
        'angular_NN_NDelta': 'Detailed balance for reactions where non-isotropic angular distributions are implemented. Angular distributions in forward and reverse reactions are not consistent.',
        'KN_KDelta': '-',
        'N_1440': 'Huge N&#963 cross-sections are problematic. Increasing the maximum cross-section to 2000 mb might improve it. See also issue #3206.',
        'N_pi_Delta_12': '-',
        'N_pi_Delta_22': '-',
        'N_pi_Delta_all': 'There seems to be a physical reason for deviations. Forbidding a special kind of spurious collisions (see issue #3025) or increasing Ntest cures it.',
        'pi_rho': '-',
        'pi_rho_f2': '-',
        'pi_rho_omega': '-',
        'pi_rho_sigma': '-',
        'pi_sigma': '-',
        'thermal_box': 'Shows only most violating reactions. Huge cross-sections are problematic. Increasing the maximum cross-section will partially improve it.',
        'Strangeness': 'Shows only most violating reactions.',
        'N_pi_deutron': '-',
        'multi_pi2_eta_etaprime': '-',
        'multi_pi3_omega': '-',
        'multi_pi3_phi': '-',
        'multi_N_pi_deutron': '-',
    }, 'angular_distributions': {
        'np_1.194': '-',
        'pp_1.25': '-',
        'pp_1.75': '-',
        'pp_2.80': '-',
    }, 'dileptons': {
        'ArKCl_1.76_filtered': 'SMASH neglects medium modifications, deviations from experimental data are expected. See <a href="https://arxiv.org/pdf/1711.10297.pdf"> dilepton paper</a> for details.',
        'CC_2.0_filtered': 'Too few dileptons are produced in np collisions which results in an underestimation in the intermedite energy region. See <a href="https://arxiv.org/pdf/1711.10297.pdf"> dilepton paper</a> for details.',
        'pp_1.25_filtered': '-',
        'pp_3.5_filtered': '-',
        'pNb_3.5_filtered': '-',
        'ArKCl_1.76_unfiltered': '-',
        'CC_2.0_unfiltered': '-',
        'pp_1.25_unfiltered': '-',
        'pp_3.5_unfiltered': '-',
        'pNb_3.5_unfiltered': '-',
    }, 'densities': {
        'AuAu_central_cell': '-',
        'pi_k_n_lambda': '-',
        'pi_rho': '-'
    }, 'energy_scan': {
        'total_multiplicity': '-',
        'midrapidity_yield':'-',
        'meanmt': '-',
        'meanpt': 'The cross section scaling factors increase linearly below &radic;<span style="text-decoration: overline">s</span><sub>NN </sub> = 200 Gev and follow a step function above &radic;<span style="text-decoration: overline">s</span><sub>NN </sub> = 200 GeV to enable running at high beam energies. This explains the abrpubt change at &radic;<span style="text-decoration: overline">s</span><sub>NN </sub> = 200 GeV.',
        'mtspectra': '-',
        'yspectra': '-',
    }, 'total_multiplicity': {'pi0' : '-', 'pi' : 'Pion multiplicities at low energies known to deviate.\
        For &#960<sup>+</sup> in pp caused by the overshoot in the inclusive pp cross section for &#960<sup>+</sup> production (see cross section pp plots).',
        'kaon' : '-', 'proton' : '-', 'lambda' : '-', 'xi' : '-', 'omega': '-'},
        'midrapidity_yield': {'pi0' : '-', 'pi' : '-', 'kaon' : '-', 'proton' : '-', 'lambda' : '-', 'xi' : '-', 'omega': '-'},
        'meanmt': {'pi0' : '-', 'pi' : '-', 'kaon' : '-', 'proton' : '-', 'lambda' : '-', 'xi' : '-', 'omega': '-'},
        'meanpt': {'pi0' : '-', 'pi' : '-', 'kaon' : '-', 'proton' : '-', 'lambda' : '-', 'xi' : '-', 'omega': '-'},
        'mtspectra': {'pi0' : '-', 'pi' : '-', 'kaon' : '-', 'proton' : '-', 'lambda' : '-', 'xi' : '-', 'omega': '-'},
        'yspectra': {'pi0' : '-', 'pi' : '-', 'kaon' : '-', 'proton' : '-', 'lambda' : '-', 'xi' : '-', 'omega': '-'},
}

descriptions_targetpage = {
   'cross_sections': 'Collection of cross sections for different elementary scattering processes. Inclusive and exclusive cross \
                      sections are presented.',
   'detailed_balance': 'Detailed balance tests in an infinite matter simulation (box) with different particle species. <br> \
                        On the left, the evolution of the box towards chemical equilibrium is displayed for all particle species considered. <br> \
                        The upper plot on the right shows the number of reactions binned in invariant mass and normalized by the length \
                        of the timestep (endtime - onset of equilibrium). <br> \
                        The lower plot on the right shows the numbers of most violating forward (tip to the right) and backward (tip to the left) reactions \
                        encountered in the corresponding setups. If detailed balance was perfectly conserved, the inwards and outwards arrows of each reaction \
                        should lie exactly on top of each other. To increase readabiliy, the total number of reactions for each channel is normalized by the average \
                        reaction number for the considered isospin group. Therefore, the normalized values are expected to lie at 1.0.',
   'energy_scan': 'Multiplicities, midrapidity yields, transverse momentum and rapidity spectra for the most abundant hadron species over a large range of beam energies \
                   from &radic;<span style="text-decoration: overline">s</span><sub>NN </sub> = 2.695 - 7000 GeV.',
   'elastic_box': 'Analysis of scattering rates in infinite matter simulations with only elastic collisions. \
                   Sensitivity of scattering rates to variation of different setup parameters is investigated. <br> \
                   The number of collisions in SMASH is, in the following plots, normalized to the theoretically expected number of collisions of a dilute gas: \
                   <br> N<sub>events</sub> t &#963  N<sup>2</sup>  N<sub>test</sub>  &lt;v<sub>rel</sub>&gt; / (2 V) <br> \
                   We therefore expect this value at &#8776; 1.0.',
   'FOPI_pions': 'Pion production compared to experimental data by the FOPI collaboration. See <a href="http://inspirehep.net/record/1471536"> SMASH paper</a> (Sections IV) for details. \
                   Potentials, Fermi motion, Pauli blocking and test particles are employed in the SMASH simulation.',
   'angular_distributions': 'Non-isotropic angular distributions for NN &#10231 N&#916 scattering. <br> For details on the implementation, \
                             see the <a href="http://inspirehep.net/record/1471536"> SMASH paper</a> (Sections II D 8 and III B).',
   'dileptons': 'Invariant mass spectra for dileptons from different elementary and heavy-ion collisions compared to HADES data.',
   'densities' : 'Time evolution of the baryon density in the central cell of a Au + Au collision and of the particle density in infinite matter simulations (box) for different degrees of freedom.',
   'total_multiplicity': 'Total particle multiplicities for &#960, K, P, &#923, &#926 and &#937 multiplets as a function of center-of-mass energy in pp and AuAu/PbPb collisions. <br> <br> \
                         Note: Kink in the region of &radic;<span style="text-decoration: overline">s</span><sub>NN </sub> &#x2248 3 - 5 GeV is due to uncertainties in the transition from resonances to strings. <br> \
                         Note: &#931<sup>0</sup> contribution is added to &#923 in AuAu and PbPb collisions.',
   'midrapidity_yield':'Particle yields at midrapidity for &#960, K, P, &#923, &#926 and &#937 multiplets as a function of center-of-mass energy in pp and AuAu/PbPb collisions. <br> <br> \
                         Note: Kink in the region of &radic;<span style="text-decoration: overline">s</span><sub>NN </sub> &#x2248 3 - 5 GeV is due to uncertainties in the transition from resonances to strings. <br> \
                         Note: &#931<sup>0</sup> contribution is added to &#923 in AuAu and PbPb collisions.',
   'meanmt': 'Mean transverse mass at midrapidity for &#960, K, P, &#923, &#926 and &#937 multiplets as a function of center-of-mass energy in pp and AuAu/PbPb collisions. <br> <br> \
                         Note: Peak in the region of &radic;<span style="text-decoration: overline">s</span><sub>NN </sub> &#x2248 3 - 5 GeV is due to uncertainties in the transition from resonances to strings.',
   'meanpt': 'Mean transverse momentum at midrapidity for &#960, K, P, &#923, &#926 and &#937 multiplets as a function of center-of-mass energy in pp and AuAu/PbPb collisions. <br> <br> \
                         Note: Peak in the region of &radic;<span style="text-decoration: overline">s</span><sub>NN </sub> &#x2248 3 - 5 GeV is due to uncertainties in the transition from resonances to strings.',
   'mtspectra': 'Transverse mass spectra for &#960, K, P, &#923, &#926 and &#937 multiplets at AGS, SPS and RHIC/LHC energies. <br> <br> \
   Note: &#931<sup>0</sup> contribution is added to &#923 in AuAu and PbPb collisions.',
   'yspectra': 'Rapidity spectra for &#960, K, P, &#923, &#926 and &#937 multiplets at AGS, SPS and RHIC/LHC energies in pp and AuAu/PbPb collisions. <br> <br> \
   Note: &#931<sup>0</sup> contribution is added to &#923 in AuAu and PbPb collisions.',
}

descriptions_frontpage = {
   'cross_sections': 'Collection of cross sections for different elementary scattering processes. Inclusive and exclusive cross \
                      sections are presented.',
   'detailed_balance': 'Detailed balance tests in an infinite matter simulation (box) with different particle species. The evolution towards chemical equilibrium \
                        and the reaction numbers of forward and backward reactions are displayed.  ',
   'energy_scan': 'Multiplicities, midrapidity yields, transverse momentum and rapidity spectra for the most abundant hadron species over a large range of beam energies from &radic;<span style="text-decoration: overline">s</span><sub>NN </sub> = 2.695 - 7000 GeV.',
   'elastic_box': 'Analysis of scattering rates in infinite matter simulations with only elastic collisions. \
                   Sensitivity of scattering rates to variation of different setup parameters is investigated.',
   'FOPI_pions': 'Pion production compared to experimental data by the FOPI collaboration.',
   'angular_distributions': 'Non-isotropic angular distributions for NN &#10231 N&#916 scattering.',
   'dileptons': 'Invariant mass spectra for dileptons from different elementary and heavy-ion collisions compared to HADES data.',
   'densities': 'Time evolution of the densities in collider and infinite matter simulations.'

}

target_headers = {'cross_sections': {
        'kplus_neutron': 'K<sup>+</sup> + n',
        'kplus_proton': 'K<sup>+</sup> + p',
        'kminus_neutron': 'K<sup>-</sup> + n',
        'kminus_proton': 'K<sup>-</sup> + p',
        'neutron_proton': 'n + p',
        'proton_proton': 'p + p',
        'antiproton_proton': 'p + p&#772',
        'piplus_proton': '&#960<sup>+</sup> + p',
        'piplus_piminus': '&#960<sup>+</sup> + &#960<sup>-</sup>',
        'piminus_proton': '&#960<sup>-</sup> + p',
    }, 'elastic_box': {
        'scatrate_vs_N': 'Particle Number Dependence',
        'scatrate_vs_Ntest': 'Test Particles Dependence',
        'scatrate_vs_T': 'Temperature Dependence',
        'scatrate_vs_V': 'Volume Dependence',
        'scatrate_vs_dt': 'Time Step Size Dependence',
        'scatrate_vs_sigma': 'Elastic Cross Section Dependence',
    }, 'FOPI_pions': '-',
       'detailed_balance': {
        'angular_NN_NDelta': 'Non-isotropically Distributed: N + N &#10231 &#916 + N, &emsp; N + N &#10231 &#916 + &#916',
        'KN_KDelta': 'K + N  &#10231 K + &#916',
        'N_1440': 'N + N(1440) &#10231 N + N, &emsp; N(1440) &#10231 &#960 + N, &emsp; N(1440) &#10231 &#960 + N',
        'N_pi_Delta_12': 'N + &#960 &#10231 &#916',
        'N_pi_Delta_22': 'N + N &#10231 &#916 + N, &emsp; N + N &#10231 &#916 + &#916 ',
        'N_pi_Delta_all': 'N + &#960 &#10231 &#916, &emsp; N + N &#10231 &#916 + N, &emsp; N + N &#10231 &#916 + &#916',
        'pi_rho': '&#961 &#10231 &#960 + &#960',
        'pi_rho_f2': '&#961 + &#961 &#10231 f<sub>2</sub>, &emsp; &#961 &#10231 &#960 + &#960, &emsp; &#960 + &#960 &#10231 f<sub>2</sub> ',
        'pi_rho_omega': '&#969 &#10231 &#960 + &#961',
        'pi_rho_sigma': '&#961 &#10231 &#960 + &#960, &emsp; &#963 &#10231 &#960 + &#960',
        'pi_sigma': '&#963 &#10231 &#960 + &#960',
        'thermal_box': 'Full Thermal Box with All Mesons and Baryons',
        'Strangeness': 'Strangeness Exchange',
        'N_pi_deutron': 'Deuteron: N + N &#10231 d&#8242;, &emsp; N + d &#10231 N + d&#8242, &emsp; &#960 + d &#10231 N + N, &emsp; &#960 + d &#10231 &#960 + d&#8242' ,
        'multi_pi2_eta_etaprime': 'Multiparticle: &#951&#8242 &#10231 &#960 + &#960 + &#951 ' ,
        'multi_pi3_omega': 'Multiparticle: &#969 &#10231 &#960 + &#960 + &#960' ,
        'multi_pi3_phi': 'Multiparticle: &#966 &#10231 &#960 + &#960 + &#960' ,
        'multi_N_pi_deutron': 'Multiparticle: d &#10231 N + &#960' ,
    }, 'angular_distributions': {
        'np_1.194': 'n + p @ p<sub>lab</sub> = 1.194 GeV',
        'pp_1.25': 'p + p @ p<sub>lab</sub> = 1.25 GeV',
        'pp_1.75': 'p + p @ p<sub>lab</sub> = 1.75 GeV',
        'pp_2.80': 'p + p @ p<sub>lab</sub> = 2.80 GeV',
    }, 'dileptons': {
        'ArKCl_1.76_filtered': 'Ar + KCl @ E<sub>kin</sub> = 1.76 GeV',
        'CC_2.0_filtered': 'C + C @ E<sub>kin</sub> = 2.0 GeV',
        'pp_1.25_filtered': 'p + p @ E<sub>kin</sub> = 1.25 GeV',
        'pp_3.5_filtered': 'p + p @ E<sub>kin</sub> = 3.5 GeV',
        'pNb_3.5_filtered': 'p + Nb @ E<sub>kin</sub> = 3.5 GeV',
        'ArKCl_1.76_unfiltered': 'Ar + KCl @ E<sub>kin</sub> = 1.76 GeV',
        'CC_2.0_unfiltered': 'C + C @ E<sub>kin</sub> = 2.0 GeV',
        'pp_1.25_unfiltered': 'p + p @ E<sub>kin</sub> = 1.25 GeV',
        'pp_3.5_unfiltered': 'p + p @ E<sub>kin</sub> = 3.5 GeV',
        'pNb_3.5_unfiltered': 'p + Nb @ E<sub>kin</sub> = 3.5 GeV',
    }, 'densities' : {
        'AuAu_central_cell': 'Au + Au @ E<sub>kin</sub> = 1.23 GeV',
        'pi_k_n_lambda': 'Thermal box at T = 160 MeV with &#961, K, N, &#916',
        'pi_rho': 'Thermal box at T = 160 MeV with &#960, &#961'
    }, 'energy_scan': {
        'total_multiplicity': 'Total Multiplicity',
        'midrapidity_yield':'Midrapidity Yield',
        'meanmt': 'Mean m<sub>T</sub>',
        'meanpt': 'Mean p<sub>T</sub>',
        'mtspectra': 'm<sub>T</sub> Spectra',
        'yspectra': 'y Spectra',
    }, 'total_multiplicity':
        {'pi0' : '&#960<sup>0</sup>',
         'pi' : '&#960<sup>&#177</sup>',
         'kaon' : 'K<sup>&#177</sup>',
         'proton' : 'p/p&#773',
         'lambda' : ' &#923/&#923&#773',
         'xi' : ' &#926<sup>&#177</sup>',
         'omega': '&#937<sup>&#177</sup>'
    }, 'midrapidity_yield':
        {'pi0' : '&#960<sup>0</sup>',
         'pi' : '&#960<sup>&#177</sup>',
         'kaon' : 'K<sup>&#177</sup>',
         'proton' : 'p/p&#773',
         'lambda' : ' &#923/&#923&#773',
         'xi' : ' &#926<sup>&#177</sup>',
         'omega': '&#937<sup>&#177</sup>'
    }, 'meanmt':
        {'pi0' : '&#960<sup>0</sup>',
         'pi' : '&#960<sup>&#177</sup>',
         'kaon' : 'K<sup>&#177</sup>',
         'proton' : 'p/p&#773',
         'lambda' : ' &#923/&#923&#773',
         'xi' : ' &#926<sup>&#177</sup>',
         'omega': '&#937<sup>&#177</sup>'
    }, 'meanpt':
        {'pi0' : '&#960<sup>0</sup>',
         'pi' : '&#960<sup>&#177</sup>',
         'kaon' : 'K<sup>&#177</sup>',
         'proton' : 'p/p&#773',
         'lambda' : ' &#923/&#923&#773',
         'xi' : ' &#926<sup>&#177</sup>',
         'omega': '&#937<sup>&#177</sup>'
    }, 'mtspectra':
        {'pi0' : '&#960<sup>0</sup>',
         'pi' : '&#960<sup>&#177</sup>',
         'kaon' : 'K<sup>&#177</sup>',
         'proton' : 'p/p&#773',
         'lambda' : ' &#923/&#923&#773',
         'xi' : ' &#926<sup>&#177</sup>',
         'omega': '&#937<sup>&#177</sup>'
    }, 'yspectra':
        {'pi0' : '&#960<sup>0</sup>',
         'pi' : '&#960<sup>&#177</sup>',
         'kaon' : 'K<sup>&#177</sup>',
         'proton' : 'p/p&#773',
         'lambda' : ' &#923/&#923&#773',
         'xi' : ' &#926<sup>&#177</sup>',
         'omega': '&#937<sup>&#177</sup>'
        },
}

energy_scan_sorted_observables = ['total_multiplicity', 'midrapidity_yield', 'meanmt', 'meanpt', 'mtspectra', 'yspectra']

E_scan_observable_dict = {'total_multiplicity' : 'Total Multiplicity',
                          'midrapidity_yield' : 'Midrapidity Yield',
                          'meanmt' : 'Mean m<sub>T</sub>',
                          'meanpt' : 'Mean p<sub>T</sub>',
                          'mtspectra' : 'm<sub>T</sub> Spectra',
                          'yspectra' : 'Rapidity Spectra',
}

energy_scan_sorted_spectra = ['pp_AGS', 'AuAuPbPb_AGS', 'pp_SPS', 'AuAuPbPb_SPS', 'pp_RHIC_LHC', 'AuAuPbPb_RHIC_LHC']

sorted_subtargets = {'angular_distributions' : ['np_1.194', 'pp_1.25', 'pp_1.75', 'pp_2.80'],
                     'detailed_balance' : ['thermal_box', 'pi_rho', 'pi_rho_f2', 'pi_rho_omega', 'pi_rho_sigma', 'pi_sigma',\
                                        'N_pi_Delta_12', 'N_pi_Delta_22', 'N_pi_Delta_all', 'N_1440', 'KN_KDelta', 'angular_NN_NDelta', \
                                        'Strangeness', 'N_pi_deutron', 'multi_pi2_eta_etaprime', 'multi_pi3_omega', 'multi_pi3_phi', 'multi_N_pi_deutron'],
                     'dileptons' : ['pp_1.25_filtered', 'pp_1.25_unfiltered', 'pp_3.5_filtered', 'pp_3.5_unfiltered', 'pNb_3.5_filtered', 'pNb_3.5_unfiltered', 'CC_2.0_filtered', 'CC_2.0_unfiltered', 'ArKCl_1.76_filtered', 'ArKCl_1.76_unfiltered'],
                     'cross_sections' : ['proton_proton', 'neutron_proton', 'antiproton_proton', 'piplus_proton', 'piminus_proton', 'piplus_piminus', 'kplus_proton', 'kminus_proton', 'kplus_neutron', 'kminus_proton'],
                     'elastic_box' : ['scatrate_vs_N', 'scatrate_vs_Ntest', 'scatrate_vs_T', 'scatrate_vs_V', 'scatrate_vs_dt', 'scatrate_vs_sigma'],
}

cross_sections_order = ['xs_process_type', 'xs_grouped', 'xs_final_grouped', 'xs_individual']

pdgs_sorted = ['111', '211', '321', '2212', '3122', '3312', '3334']

pdgs_to_name = {'111' : 'pi0', '211' : 'pi', '321': 'kaon', '2212': 'proton',
                '3122': 'lambda', '3312': 'xi', '3334': 'omega'}
