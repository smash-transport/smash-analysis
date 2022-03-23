import sys
import numpy as np
import argparse
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mult_and_spectra import BulkObservables
from collections import defaultdict
from itertools import cycle
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'../../python_scripts')
from txt_io import load_table
from glob import glob
import smash_basic_scripts as sb
exec(compile(open(os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts/common_plotting.py', "rb").read(), os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts/common_plotting.py', 'exec'))
import comp_to_prev_version as cpv

colours = ["darkred","deepskyblue","deeppink","limegreen","orange","red", "magenta", "chartreuse", "cyan", "darkmagenta", "yellow"]# , "lavender", "navajowhite", "midnightblue", "yellow", "blue", "gold", "darkgreen"]

def fullsplit_path(the_path):
    "Turns /a/b/c/d/e/f.dat into ['a','b','c','d','e','f.dat'] correctly for different filesystems"
    return os.path.normpath(the_path).strip(os.path.sep).split(os.path.sep)

def store_results(results_file, results_list, version, quantity):
    '''
    Store results for different quantities to allow comparisons to future
    SMASH version.
    The columns in transverse mass and rapidity spectra correspond to the
    different collision energies. For all other quantities, they correspond to
    the pdg numbers of analyzed particles.
    The generated files are stored in the same directory as the plots.
    '''

    store_file = open(results_file, "w")

    if ('spectra' not in quantity):
        title_line = 'sqrt_s' + "\t"
        title_line += '\t'.join(str(e) for e in results_list[0])

        store_file.write(title_line + '\n')
        for energy in range(0,len(results_list[1][0])):
            result_line = str(results_list[1][0][energy])
            for pdg in range(0,len(results_list[0])):
                result_line += "\t" + str(results_list[2][pdg][energy])
            store_file.write(result_line + "\n")

    else:
        if quantity == 'mtspectra':
            title_line = 'mt-m0' + "\t"
            title_line += '\t'.join(str(e) for e in results_list[0])
        elif quantity == 'yspectra':
            title_line = 'y' + "\t"
            title_line += '\t'.join(str(e) for e in results_list[0])
        elif quantity == 'ptspectra':
            title_line = 'pt' + "\t"
            title_line += '\t'.join(str(e) for e in results_list[0])
        elif quantity == 'v2spectra':
            title_line = 'pt' + "\t"
            title_line += '\t'.join(str(e) for e in results_list[0])
        else:
            print("No scheme for formatting " + str(quantity) + " found.")

        store_file.write(title_line + '\n')
        for i in range(0,len(results_list[1][0])):
            result_line = str(results_list[1][0][i])
            for energy in range(0,len(results_list[0])):
                result_line += "\t" + str(results_list[2][energy][i])
            store_file.write(result_line + "\n")

    store_file.write("# " + str(version))


def determine_collider(energy):
    if (0.0 <= energy < 5.5):
        collider = 'AGS'
    elif (5.5 <= energy < 19.6):
        collider = 'SPS'
    elif (19.6 <= energy):
        collider = 'RHIC_LHC'
    else:
        print(energy)
        raise SystemExit('Collider could not be determined.')


    return collider


class DataTree:
    """
    Implements a tree convenient to store data as a dictionary of structure:

      dict[<physical quantity>][<colliding system>][<pdg>][<energy>] = result

    """

    def __init__(self):
       self.the_dict = {}
       self.quantities = set()
       self.colliding_systems = set()
       self.pdglist = set()
       self.energies = set()
       self.midrapidity_cut = 0.5

    def add(self, the_list):
        "adding [a0, a1, a2, a3, a4, a5] to dictionary so that dict[a0][a1][a2]...[a4] = a5"
        current_leaf = self.the_dict
        assert(len(the_list) == 5)
        for x in the_list[:-2]:
            if x not in current_leaf:
                current_leaf[x] = {}
            current_leaf = current_leaf[x]
        current_leaf[the_list[-2]] = the_list[-1]
        self.quantities.add(the_list[0])
        self.colliding_systems.add(the_list[1])
        self.pdglist.add(the_list[2])
        self.energies.add(the_list[3])

    def is_in_dict(self, the_list):
        "for list [a0, a1, ..., an] checking if dict[a0][a1]...[an] exists"
        current_leaf = self.the_dict
        for x in the_list:
            if (x not in current_leaf):
                return False
            current_leaf = current_leaf[x]
        return True

    def print_summary_table(self):
        " Print the saved data in a structured manner "
        energies_sorted = sorted(list(self.energies))
        pdglist_sorted = sorted(list(self.pdglist), key = lambda x : abs(x))
        for quantity in self.quantities:
            if ('spectra' in quantity): continue
            for colliding_system in self.colliding_systems:
                if (not self.is_in_dict([quantity, colliding_system])): continue
                print(quantity, colliding_system)
                print('-' * (10*len(self.pdglist) + 18))
                print('# sqrts ', end=' ')
                for pdg in pdglist_sorted:
                    print('%10s' % pdg, end=' ')
                print()
                for energy in energies_sorted:
                    any_pdg_for_this_energy = False
                    for pdg in self.pdglist:
                        if (self.is_in_dict([quantity, colliding_system, pdg, energy])):
                            any_pdg_for_this_energy = True
                    if (not any_pdg_for_this_energy): continue
                    print('%9.4f ' % energy, end=' ')
                    for pdg in pdglist_sorted:
                        if (not self.is_in_dict([quantity, colliding_system, pdg, energy])):
                            print(' '*10, end=' ')
                            continue
                        x = self.the_dict[quantity][colliding_system][pdg][energy]
                        if isinstance(x, tuple): x = x[0]
                        print('%9.4f ' % x, end=' ')
                    print()
                print('-' * (10*len(self.pdglist) + 18))

    def read_theory(self, input_files, PbPb_as_AuAu = True):
        " Consume all the results of theoretical calculations from given directories "
        input_directories = set([os.path.dirname(f) for f in input_files])
        quantities_list = ['yspectra', 'mtspectra', 'ptspectra', 'v2spectra', 'meanmt0_midrapidity',  'meanpt_midrapidity',  'midrapidity_yield',  'total_multiplicity']
        for directory in input_directories:
            dir_split = fullsplit_path(directory)
            colliding_system = dir_split[-1]
            if (PbPb_as_AuAu and (colliding_system == 'AuAu' or colliding_system == 'PbPb')):
                colliding_system = 'AuAu/PbPb'
            try:
                energy = float(dir_split[-2])
            except Exception:
                if(dir_split[-2]=="RHIC"):
                    energy=200.0
                    colliding_system = 'afterburner'
                elif(dir_split[-2]=="LHC"):
                    energy=5000.0
                    colliding_system = 'afterburner'
                else:
                    raise SystemExit('Path could not be parsed.')
            # trunc for easier matching with experimental energies
            energy = np.trunc(energy * 1.e2) / 1.e2
            files = tuple([directory + os.path.sep + quantities_list[i] + '.txt' for i in range(len(quantities_list))])
            data = BulkObservables.read(files)
            pdglist = data.pdglist
            Npdg = len(data.pdglist)
            midrapidity_cut = data.midrapidity_cut
            for q in quantities_list:
                for i in range(Npdg):
                    # Lambdas and antilambdas need Sigma^0 to be added, because in AA collisions
                    # Lambda and Sigma^0 are indistinguishable
                    if (colliding_system != 'pp' and abs(pdglist[i]) == 3122):
                        pdg_to_find = 3212 if pdglist[i] == 3122 else -3212
                        j = pdglist.index(pdg_to_find)
                        data.total_multiplicity[i] += data.total_multiplicity[j]
                        data.midrapidity_yield[i] += data.midrapidity_yield[j]
                        data.mthist[i] += data.mthist[j]
                        data.pthist_midrapidity[i] += data.pthist_midrapidity[j]
                        data.yhist[i] += data.yhist[j]
                    if (q == 'total_multiplicity'): to_dict = float(data.total_multiplicity[i]) / data.nevents
                    elif (q == 'midrapidity_yield'): to_dict = float(data.midrapidity_yield[i]) / data.nevents
                    elif (q == 'meanpt_midrapidity'): to_dict = float(data.meanpt_midrapidity[i])
                    elif (q == 'meanmt0_midrapidity'): to_dict = float(data.meanmt0_midrapidity[i])
                    elif (q == 'mtspectra'): to_dict = (data.mtbins, \
                                                        data.mthist[i].astype(float) / data.nevents)
                    elif (q == 'yspectra'): to_dict = (data.ybins,\
                                                       data.yhist[i].astype(float) / data.nevents)
                    elif (q == 'ptspectra'): to_dict = (data.ptbins,\
                                                       data.pthist_midrapidity[i].astype(float) / data.nevents)
                    # for v2:
                    # Normalization to number of detected particles in given pT bin, to make sure
                    # events with zero entries are not considered.
                    elif (q == 'v2spectra'):
                        if args.with_v2:
                            to_dict = (data.ptbins,\
                                       data.v2[i].astype(float) / data.pthist_midrapidity[i].astype(float))
                        else: continue
                    else: print('error: unexpected quantity ', q)
                    self.add([q, colliding_system, pdglist[i], energy, to_dict])
        return data.smash_version

    def read_experiment(self, input_files_experiment):
        " Consume experimental data from the given input files. "

        quantity_exp_dict = {
          ('multiplicity', 'ally') : 'total_multiplicity',
          ('midy_yield', 'midy') : 'midrapidity_yield',
          ('mean_pT', 'midy') : 'meanpt_midrapidity',
          ('mean_mT', 'midy') : 'meanmt0_midrapidity'
        }

        colliding_system_exp_dict = {
          'AA' : 'AuAu/PbPb',
          'pp' : 'pp'
        }

        pdg_exp_dict = {
          'piplus' : 211,
          'piminus' : -211,
          'pizero' : 111,
          'kplus' : 321,
          'kminus' : -321,
          'proton' : 2212,
          'antiproton' : -2212,
          'lambda' : 3122,
          'antilambda' : -3122,
          'xi' : 3312,
          'antixi' : -3312,
          'omega' : 3334,
          'antiomega' : -3334
        }

        def exp_quantity_from_filename(folder, fname):
            for i,j in quantity_exp_dict:
                if ((i in folder) and (j in fname)):
                    return quantity_exp_dict[(i,j)]
            return None

        def pdg_from_filename(fname):
            for i in  pdg_exp_dict:
                if fname.split('_')[1] == i: return pdg_exp_dict[i]
            return None

        for exp_file in input_files_experiment:
            file_name_split = fullsplit_path(exp_file)
            expfolder, expfile = file_name_split[-3], file_name_split[-1]
            colliding_system = file_name_split[-2]
            if colliding_system not in ['pp', 'AA']: continue
            colliding_system = colliding_system_exp_dict[colliding_system]
            extension = expfile.split('.')[-1]
            expfile = ''.join(expfile.split('.')[:-1])
            if (extension != 'exp'): continue
            quantity = exp_quantity_from_filename(expfolder, expfile)
            pdg = pdg_from_filename(expfile)
            if ((quantity is None) or (colliding_system is None) or (pdg is None)):
                continue
            # print quantity, colliding_system, pdg
            with open(exp_file, 'r') as f:
                # Assuming that datafiles are small, kB, but not GB
                for line in f:
                    line_no_comment = line.split('!')[0].strip()
                    if (line_no_comment == ''): continue
                    if ('*' in line_no_comment): continue
                    l = line_no_comment.split()
                    if (len(l) == 3):
                        energy, exp_meas, exp_meas_err = [float(x) for x in l]
                    elif (len(l) == 4):
                        energy, exp_meas, err1, err2 = [float(x) for x in l]
                        exp_meas_err = np.sqrt(err1*err1 + err2*err2)
                    else:
                        print("Unexpected line: ", line)
                    self.add([quantity, colliding_system, pdg, energy, (exp_meas, exp_meas_err)])

    def read_experiment_spectra(self, input_files):
        " Read rapidity and mt spectra from given files"
        for exp_file in input_files:
            # print exp_file
            # By default it is assumed that it is rapidity spectrum, unless it has dndmt in its name
            is_dndmt = 'mt_spectrum' in exp_file
            def sqrts(elab):
                mN = 0.938
                return np.sqrt(2.0 * mN * (2.0 * mN + elab))
            for i in ['ekin', 'elab', 'elb', 'ecm']:
                if (i in exp_file):
                    energy_elab = float(exp_file.split(i)[1].split('_')[0])
                    energy = sqrts(energy_elab) if (i != 'ecm') else energy_elab
                    break
            energy = np.trunc(energy * 1.e2) / 1.e2
            name_to_pdg = {
                'piplus' : 211,
                'piminus' : -211,
                'pizero' : 111,
                'kplus' : 321,
                'kminus' : -321,
                'proton' : 2212,
                'antiproton' : -2212,
                'lambda' : 3122,
                'antilambda' : -3122,
                'xi' : 3312,
                'antixi' : -3312,
                'omega' : 3334,
                'antiomega' : -3334,
            }
            for part_name in list(name_to_pdg.keys()):
                if part_name == exp_file.split('/')[-1].split('.')[0].split('_')[-1]:
                    break

            pdg = name_to_pdg[part_name]
            # print energy, part_name, name_to_pdg[part_name]
            def is_number(s):
                if (s.strip() == ''): return False
                try:
                    float(s)
                    return True
                except ValueError:
                    return False
            nlines_to_skip = 0
            with open(exp_file, 'r') as f:
                 for line in f:
                     a = line.strip().split()
                     if (a and (all([is_number(x) for x in a]))): break
                     nlines_to_skip += 1
            data = np.loadtxt(exp_file, skiprows = nlines_to_skip, comments = ['!', '#'])
            where_to_add = 'yspectra' if not is_dndmt else 'mtspectra'
            if 'pp' in exp_file:
                self.add([where_to_add, 'pp', pdg, energy, (data[:,0], data[:,1], data[:,2])])
            else:
                self.add([where_to_add, 'AuAu/PbPb', pdg, energy, (data[:,0], data[:,1], data[:,2])])
            # print data

def plotting(data1, data2, config_file, smash_code_version, output_folder):
    plot_counter=0
    quantities = data1.quantities.union(data2.quantities)
    pdglist = data1.pdglist.union(data2.pdglist)
    pdglist_abs = np.unique(np.abs(np.array(list(pdglist))))
    colliding_systems = data1.colliding_systems.union(data2.colliding_systems)
    if('afterburner' in colliding_systems):
        colliding_systems={'afterburner'}
    energies = sorted(list(data1.energies.union(data2.energies)))
    for quantity in quantities:
        if ('spectra' in quantity): continue
        #if ('afterburner' in colliding_systems): continue
        collected_results_pp = [[],[],[]]
        collected_results_AuAuPbPb = [[],[],[]]
        collected_results_afterburner = [[],[],[]]
        for pdg_abs in pdglist_abs:
            # Sigma0 is added to Lambda plots
            if (pdg_abs == 3212): continue
            # We do not want the Deuteron plots to be displayed, because SMASH
            # needs to be modified for useful results (cross section cut off)
            if (pdg_abs == 1000010020): continue
            pdg_one_sort = []
            if (pdg_abs in pdglist): pdg_one_sort.append(pdg_abs)
            if (-pdg_abs in pdglist): pdg_one_sort.append(-pdg_abs)
            for pdg in pdg_one_sort:
                for colliding_system in colliding_systems:
                    if (pdg == pdg_abs):
                        plot_format = '-'
                        plot_label = colliding_system
                        exp_fmt = 'o'
                    else:
                        plot_format = '--'
                        plot_label = ''
                        exp_fmt = 's'

                    if (colliding_system == 'pp'):
                        linewidth = 5
                        plot_color = 'midnightblue'
                    else:
                        linewidth = 5
                        plot_color = 'darkred'
                    if data1.is_in_dict([quantity, colliding_system, pdg]):
                        theory_vs_energy = data1.the_dict[quantity][colliding_system][pdg]
                        x, y  = list(zip(*sorted(theory_vs_energy.items())))
                        plt.plot(x, y, plot_format, label = plot_label, color = plot_color, lw = linewidth)

                        # store results in list to later save it for future comparison
                        if (colliding_system == 'pp'):
                            collected_results_pp[0].append(pdg)
                            if (collected_results_pp[1] == []): collected_results_pp[1].append(x)
                            collected_results_pp[2].append(y)
                            if args.comp_prev_version:
                                import comp_to_prev_version as cpv
                                # read and plot results from previous version
                                filename_prev = quantity + '_' + colliding_system.replace('/', '')
                                prev_SMASH_version = cpv.plot_previous_results('energy_scan', '', filename_prev + '.txt',
                                                    plot_color = plot_color, pdg = pdg, plot_style = plot_format)
                        if (colliding_system == 'AuAu/PbPb'):
                            collected_results_AuAuPbPb[0].append(pdg)
                            if (collected_results_AuAuPbPb[1] == []): collected_results_AuAuPbPb[1].append(x)
                            collected_results_AuAuPbPb[2].append(y)
                            if args.comp_prev_version:
                                import comp_to_prev_version as cpv
                                # read and plot results from previous version
                                filename_prev = quantity + '_' + colliding_system.replace('/', '')
                                prev_SMASH_version = cpv.plot_previous_results('energy_scan', '', filename_prev + '.txt',
                                                    plot_color = plot_color, pdg = pdg, plot_style = plot_format)
                        if (colliding_system == 'afterburner'):
                            collected_results_afterburner[0].append(pdg)
                            if (collected_results_afterburner[1] == []): collected_results_afterburner[1].append(x)
                            collected_results_afterburner[2].append(y)
                            if args.comp_prev_version:
                                import comp_to_prev_version as cpv
                                # read and plot results from previous version
                                filename_prev = quantity  +  '_' + colliding_system.replace('/', '')
                                prev_SMASH_version = cpv.plot_previous_results('afterburner', '', filename_prev + '.txt',
                                                    plot_color = plot_color, pdg = pdg, plot_style = plot_format)

                    if (colliding_system != 'afterburner' and data2.is_in_dict([quantity, colliding_system, pdg])):
                        exp_vs_energy = data2.the_dict[quantity][colliding_system][pdg]
                        x_exp, y_exp = list(zip(*sorted(exp_vs_energy.items())))
                        y_exp_values, y_exp_errors = list(zip(*y_exp))
                        plt.errorbar(x_exp, y_exp_values, yerr = y_exp_errors, fmt = exp_fmt,\
                                     color = plot_color, elinewidth = 2, markeredgewidth = 0)
                    elif (colliding_system == 'afterburner' and data2.is_in_dict([quantity, 'AuAu/PbPb', pdg])):
                        exp_vs_energy = data2.the_dict[quantity]['AuAu/PbPb'][pdg]
                        x_exp, y_exp = list(zip(*sorted(exp_vs_energy.items())))
                        y_exp_values, y_exp_errors = list(zip(*y_exp))
                        plt.errorbar(x_exp, y_exp_values, yerr = y_exp_errors, fmt = exp_fmt,\
                                     color = plot_color, elinewidth = 2, markeredgewidth = 0)


            if args.comp_prev_version:
                #dummy, for combined legend entry of previous results.
                plt.plot(1, 0.0, linestyle = '-', linewidth = 10, zorder = 1,
                        alpha = 0.2, color='dimgrey',label=prev_SMASH_version)

            plt.xlabel('$\sqrt{s_{NN}} [GeV]$')
            plt.xscale('log')

            if (quantity in ['total_multiplicity', 'midrapidity_yield']):
                if np.all(y == 0):
                    print('No positive values encountered in ' + str(quantity) + ' for ' + str(pdg) +\
                          '. Cannot log-scale the y axis, scale will be linear.')
                else:
                     plt.yscale('log', nonposy='clip')
            if( quantity in ['total_multiplicity', 'midrapidity_yield', 'meanmt0_midrapidity', 'meanpt_midrapidity']  and 'afterburner' in colliding_systems):
                plt.xlim([190,5100])
                plt.text(0.5, 0.5, 'SMASH-vHLLE-hybrid', fontsize=40, color='gray', ha='right', va='bottom', alpha=0.5, transform=plt.gca().transAxes)
            hadron_name = sb.pdg_to_name(pdg_abs, config_file)
            antihadron_name = sb.pdg_to_name(-pdg_abs, config_file)
            plot_title = hadron_name
            if not pdg_abs in [111, 333]:  # antiparticle of itself
                plot_title += ' (' + antihadron_name + ' dashed, squares)'
            title_dict = {
                'total_multiplicity' : ' $4\pi$ multiplicity',
                'midrapidity_yield' :  ' $dN/dy|_{y = 0}$',
                'meanmt0_midrapidity' : ' $<m_{T}>|_{y = 0}$, $m_{T} = \sqrt{p_T^2 + m^2} - m$',
                'meanpt_midrapidity' :  ' $<p_{T}>|_{y = 0}$'
            }
            plt.title(plot_title)
            plt.ylabel(title_dict[quantity])
            plt.legend()
            plt.figtext(0.8, 0.94, " SMASH code:      %s\n SMASH analysis: %s" % \
                         (smash_code_version, sb.analysis_version_string()), \
                         color = "gray", fontsize = 10)
            plt.savefig(output_folder + '/' + quantity + str(pdg_abs) + '.pdf')
            plt.clf()

        # Save results plotted above for future comparison
        if('afterburner' in colliding_systems):
            filename_afterburner = quantity + '_' + 'afterburner' + '.txt'
            store_results(output_folder + '/' + filename_afterburner, collected_results_afterburner, smash_code_version, quantity)
        else:
            filename_AuAuPbPb = quantity + '_' + 'AuAuPbPb' + '.txt'
            filename_pp = quantity + '_' + 'pp' + '.txt'
            store_results(output_folder + '/' + filename_AuAuPbPb, collected_results_AuAuPbPb, smash_code_version, quantity)
            store_results(output_folder + '/' + filename_pp, collected_results_pp, smash_code_version, quantity)
    # Plotting spectra, only those, where some data is present
    for quantity in quantities:
        if (quantity not in ['yspectra', 'mtspectra', 'ptspectra', 'v2spectra']): continue
        if (quantity == 'v2spectra' and not args.with_v2): continue
        for pdg in pdglist:
            if (abs(pdg) == 3212): continue
            if (abs(pdg) == 1000010020): continue
            collected_results_pp = [[],[],[]]
            collected_results_AuAuPbPb = [[],[],[]]
            collected_results_afterburner=[[],[],[]]
            for colliding_system in colliding_systems:
                #if not data2.is_in_dict([quantity, colliding_system, pdg]): continue
                #if not data1.is_in_dict([quantity, colliding_system, pdg]): continue

                # colors list for plotting
                col = cycle(colours)
                # to scale curves in mT and pT spectra by powers of 10 -> readability
                scaling_counter = -1
                for element, energy in enumerate(energies):
                    collider = determine_collider(energy)
                    in_theory = data1.is_in_dict([quantity, colliding_system, pdg, energy])
                    if colliding_system=='afterburner':
                        in_experiment = data2.is_in_dict([quantity, 'AuAu/PbPb', pdg, energy])
                        if not in_theory:
                            in_experiment=False
                            continue

                    else:
                        in_experiment = data2.is_in_dict([quantity, colliding_system, pdg, energy])
                    if (in_experiment and not in_theory):
                        print(energy, colliding_system, pdg, in_theory, in_experiment, \
                              ': there is experimental data, but no SMASH calculation!')
                    if( not in_experiment and not in_theory):
                        continue
                    if (in_theory):
                        plot_color = next(col)
                        bin_edges, y = data1.the_dict[quantity][colliding_system][pdg][energy]
                        bin_edges = np.array(bin_edges)
                        bin_width = bin_edges[1:] - bin_edges[:-1]
                        x = 0.5*(bin_edges[1:] + bin_edges[:-1])
                        y = np.array(y)
                        # dN/dy
                        if (quantity == 'yspectra'):
                            y /= bin_width
                            plt.plot(x, y, '-', lw = 4, color = plot_color, label = str(energy))
                        # dN/dmT
                        pole_masses = {111  : 0.138,
                                       211  : 0.138,
                                       321  : 0.495,
                                       2212 : 0.938,
                                       3122 : 1.116,
                                       3312 : 1.321,
                                       3334 : 1.672,
                                       1000010020 : 1.8756,
                                       3212 : 1.189 }
                        m0 = pole_masses[abs(pdg)]
                        if (quantity == 'mtspectra'):
                            if(colliding_system!='afterburner'):
                                scaling_counter += 1
                            else:
                                scaling_counter=0
                            y /= ((x + m0) * bin_width) * (2.0 * data1.midrapidity_cut)  # factor 2 because [-y_cut; y_cut]
                            if np.all(y == 0):          # rescale y-axis to be linear if mtspectra of current energy are 0, but those
                                plt.yscale('linear')    # of the previous energy were not, so that the scale was already set to log scale.
                            plt.plot(x, y * 10**scaling_counter, '-', lw = 4, color = plot_color,
                                    label = str(energy) + r' $\times \ $10$^{\mathsf{' + str(scaling_counter) + r'}}$')
                        # dN/dpT
                        if (quantity == 'ptspectra'):
                            if(colliding_system!='afterburner'):
                                scaling_counter += 1
                            else:
                                scaling_counter=0
                            y /= (bin_width * x) * (2.0 * data1.midrapidity_cut)  # factor 2 because [-y_cut; y_cut]
                            if np.all(y == 0):          # rescale y-axis to be linear if ptspectra of current energy are 0, but those
                                plt.yscale('linear')    # of the previous energy were not, so that the scale was already set to log scale.
                            plt.plot(x, y * 10**scaling_counter, '-', lw = 4, color = plot_color,
                                    label = str(energy) + r' $\times \ $10$^{\mathsf{' + str(scaling_counter) + r'}}$')
                        # v2
                        if (quantity == 'v2spectra'):
                            y /= (bin_width) * (2.0 * data1.midrapidity_cut)  # factor 2 because [-y_cut; y_cut]
                            plt.plot(x, y, '-', lw = 4, color = plot_color,
                                    label = str(energy))


                        # store results in list to later save for future comparison
                        if (colliding_system == 'pp'):
                            collected_results_pp[0].append(energy)
                            if (collected_results_pp[1] == []):collected_results_pp[1].append(x)
                            collected_results_pp[2].append(y)
                        if (colliding_system == 'AuAu/PbPb'):
                            collected_results_AuAuPbPb[0].append(energy)
                            if (collected_results_AuAuPbPb[1] == []): collected_results_AuAuPbPb[1].append(x)
                            collected_results_AuAuPbPb[2].append(y)
                        if (colliding_system == 'afterburner'):
                            collected_results_afterburner[0].append(energy)
                            if (collected_results_afterburner[1] == []): collected_results_afterburner[1].append(x)
                            collected_results_afterburner[2].append(y)

                        # read and plot results from previous version
                        if args.comp_prev_version and quantity != 'v2spectra':
                            import comp_to_prev_version as cpv
                            #v2 is not regularly run, old results are neither produced nor stored

                            if(colliding_system == 'afterburner'):
                                filename_prev = quantity + '_' + colliding_system.replace('/', '') + str(pdg)
                                prev_SMASH_version =  cpv.plot_previous_results('afterburner', '', filename_prev + '.txt',
                                                  energy = energy, plot_color = plot_color, scaling_counter = 0)
                            else:
                                filename_prev = quantity + '_' + colliding_system.replace('/', '') + '_' + str(pdg)
                                prev_SMASH_version =  cpv.plot_previous_results('energy_scan', '', filename_prev + '.txt',
                                                  energy = energy, plot_color = plot_color, scaling_counter = scaling_counter)

                    if (in_experiment):
                        if (colliding_system != 'afterburner'):
                            x, y, y_err = data2.the_dict[quantity][colliding_system][pdg][energy]
                        else:
                            x, y, y_err = data2.the_dict[quantity]['AuAu/PbPb'][pdg][energy]
                            scaling_counter=0
                        if (quantity == 'mtspectra'):
                            plt.errorbar(x, y * 10**scaling_counter, yerr = y_err, fmt = 'o', color = plot_color)
                        elif (quantity == 'ptspectra'):
                            plt.errorbar(x, y * 10**scaling_counter, yerr = y_err, fmt = 'o', color = plot_color)
                        elif (quantity == 'v2spectra'):
                            plt.errorbar(x, y, yerr = y_err, fmt = 'o', color = plot_color)
                        else: # yspectra
                            plt.errorbar(x, y, yerr = y_err, fmt = 'o', color = plot_color)
                    title_dict = {
                        'yspectra' : '$dN/dy$',
                        'mtspectra' : '$1/m_{T} \ d^2N/dm_{T} dy$ [GeV$^{-2}$]',
                        'ptspectra' : '$1/p_{T} \ d^2N/dp_{T} dy$ [GeV$^{-2}$]',
                        'v2spectra' : '$v_2$',
                    }

                    plot_title = sb.pdg_to_name(pdg, config_file)
                    plot_title += ' in ' + colliding_system + ' collisions'
                    if colliding_system=='afterburner' :
                        plot_title =  sb.pdg_to_name(pdg, config_file)+' in AuAu/PbPb collisions'
                    plt.title(plot_title)
                    plt.figtext(0.15, 0.94, " SMASH code:      %s\n SMASH analysis: %s" % \
                        (smash_code_version, sb.analysis_version_string()), \
                        color = "gray", fontsize = 10)
                    if (quantity == 'mtspectra' or quantity == 'ptspectra'):
                        if np.all(y == 0):
                            print('No positive values encountered in ' + str(quantity) + ' for ' + str(pdg) +\
                                  '. Cannot log-scale the y axis, scale will be linear.')
                        else:
                             plt.yscale('log', nonposy='clip')
                        if (quantity == 'mtspectra'):
                            plt.xlabel('$m_{T} - m_{0}$ [GeV]')
                        else:
                            plt.xlabel('$p_{T}$ [GeV]')
                    elif (quantity == 'yspectra'):
                        plt.xlabel('$y$')
                    else:
                        plt.xlabel('$p_{T}$ [GeV]')

                    plt.ylabel(title_dict[quantity])
                    plot_counter=plot_counter+1
                    if (determine_collider(energy) != determine_collider(energies[(element + 1) % len(energies)]) or (colliding_system == 'afterburner' and plot_counter==1)):
                        if args.comp_prev_version:
                            #dummy for legend entry of combined previous results.
                            import comp_to_prev_version as cpv
                            #v2 is not regularly run, old results are neither produced nor stored

                            if (colliding_system != 'afterburner'):
                                filename_prev = quantity + '_' + colliding_system.replace('/', '') + '_' + str(pdg)
                                prev_SMASH_version =  cpv.plot_previous_results('energy_scan', '', filename_prev + '.txt',
                                                  energy = energy, plot_color = 'midnightblue', scaling_counter = scaling_counter)
                            else:
                                filename_prev = quantity + '_' + colliding_system.replace('/', '') + str(pdg)
                                prev_SMASH_version =  cpv.plot_previous_results('afterburner', '', filename_prev + '.txt',
                                                  energy = energy, plot_color = 'midnightblue', scaling_counter = 0)
                                plt.text(0.5, 0.5, 'SMASH-vHLLE-hybrid', fontsize=40, color='gray', ha='right', va='bottom', alpha=0.5, transform=plt.gca().transAxes)
                            plt.plot(1,0.0, linestyle = '-', linewidth = 10, zorder = 1,
                                    color='dimgrey', label=prev_SMASH_version, alpha = 0.2)
                        plt.legend(loc= 'upper right', title = '$\sqrt{s} \ $ [GeV] =' , ncol = 1, fontsize = 26)
                        plt.savefig(output_folder + '/' + quantity + '_' + colliding_system.replace('/', '') + '_' + collider + '_' + str(pdg) + '.pdf')
                        plt.clf()
                        plt.close()
                        scaling_counter = -1   #re-initialize as generating a new plot

            # Save results plotted above for future comparison
            if('afterburner' in colliding_systems):
                filename_afterburner = quantity + '_' + 'afterburner' + str(pdg) + '.txt'
                store_results(output_folder + '/' + filename_afterburner, collected_results_afterburner, smash_code_version, quantity)
            else:
                filename_AuAuPbPb = quantity + '_' + 'AuAuPbPb'  + str(pdg)+ '.txt'
                filename_pp = quantity + '_' + 'pp' + str(pdg) + '.txt'
                store_results(output_folder + '/' + filename_AuAuPbPb, collected_results_AuAuPbPb, smash_code_version, quantity)
                store_results(output_folder + '/' + filename_pp, collected_results_pp, smash_code_version, quantity)


if __name__ == '__main__':

    """
    Typical usage:
    python mult_and_spectra_plotter.py
       --theory energy_scan/*/AuAu/*.txt energy_scan/*/PbPb/*.txt energy_scan/*/pp/*.txt
       --config_file energy_scan/2.695/AuAu/data/1/config.yaml
       --experiment ../experimental_data/multiplicity/*/*.exp
                    ../experimental_data/midy_yield/*/*.exp
                    ../experimental_data/mean_pT/*/*.exp
                    ../experimental_data/mean_mT/*/*.exp
       --experiment_spectra ../experimental_data/mt_spectrum/*/*.exp
                            ../experimental_data/y_spectrum/*/*.exp

    Printing summary of some experimental data (works for theory too):
    python mult_and_spectra_plotter.py
        --config_file energy_scan/2.695/AuAu/data/1/config.yaml
        --print_summary True
        --experiment ../experimental_data/multiplicity/*/*.exp
                    ../experimental_data/midy_yield/*/*.exp
                    ../experimental_data/mean_pT/*/*.exp
                    ../experimental_data/mean_mT/*/*.exp
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--theory", nargs = '+', required = False,
                        help = "files of structure <path>/<colliding system>/<energy_sqrts>/file.txt")
    parser.add_argument("--theory_to_compare", nargs = '+', required = False,
                        help = "folders of structure <path>/<colliding system>/<energy_sqrts>")
    parser.add_argument("--experiment", nargs = '+', required = False, default = '',
                        help = "files with experimental data")
    parser.add_argument("--experiment_spectra", nargs = '+', required = False, default = '',
                        help = "files with experimental data")
    parser.add_argument("--config_file", type = str, required = True)
    parser.add_argument("--output_folder", type = str, required = False, default = '.')
    parser.add_argument("--PbPb_as_AuAu", type = bool, required = False, default = True,
                        help = "plot AuAu and PbPb data as separate curves")
    parser.add_argument("--print_summary", type = bool, required = False, default = False,
                        help = "Print summary tables for consumed data.")
    parser.add_argument("--with_v2", type = bool, required = False, default = False,
                        help = "Also plot v2 spectra.\
                        Caution: Large statistics are required and results are \
                        only meaningful above 40 AGeV beam energies.")
    parser.add_argument("--comp_prev_version", action='store_true',
                        help = "Plot comparison to previous SMASH version.")
    args = parser.parse_args()

    theory = DataTree()
    smash_version = ''
    if (args.theory):
        smash_version = theory.read_theory(args.theory, args.PbPb_as_AuAu)
    experiment = DataTree()
    if (args.experiment):
        experiment.read_experiment(args.experiment)
    if (args.experiment_spectra):
        experiment.read_experiment_spectra(args.experiment_spectra)

    if (args.print_summary):
         theory.print_summary_table()
         experiment.print_summary_table()

    plotting(theory, experiment, args.config_file, smash_version, args.output_folder)
