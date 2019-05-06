"""
  This script plots the ratio of two spectra, which is particularly
  useful for seeing if switching on/off Pauli blocking changes distribution
  in the box from Boltzmann to Fermi-Dirac.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import smash_basic_scripts as sb

def read_spectrum(input_file):
    """ Read in the dN/dE spectra from the input_file """
    with open(input_file, 'r') as f:
        f.readline()
        total_events = int(f.readline())
        f.readline()
        tstart = float(f.readline())
        f.readline()
        pdg_list = f.readline().rstrip().split(',')
    contents = np.loadtxt(input_file, skiprows = 6)
    return {'nevents': total_events,
            't0': tstart,
            'pdg_list': pdg_list,
            'spectra': contents}

def plot_spectra_ratio(input1, input2, output_file):
    """ Writes ratio of spectra from files input1 and input2 to output_file. """
    r1 = read_spectrum(input1)
    r2 = read_spectrum(input2)
    if (not (r1['nevents'] == r2['nevents'] and
        r1['t0'] == r2['t0'] and
        r1['pdg_list'] == r2['pdg_list'])):
        sys.exit("Files have to be obtained for equal conditions")

    plt.xlabel("E, GeV", fontsize = '20')
    plt.ylabel("$\\frac{dN_{Pauli}/dE}{dN_{no Pauli}/dE}$", fontsize = '20')
    plt.title("Pauli blocking: spectra ratio, %d events" % r1['nevents'], fontsize = '20')

    for i in xrange(len(r1['pdg_list'])):
        plt.errorbar(r1['spectra'][2*i], r1['spectra'][2*i+1]/r2['spectra'][2*i+1],
                     yerr = np.sqrt(2.0/r1['spectra'][2*i+1]),
                     label=sb.pdg_to_name(int(r1['pdg_list'][i]))+", m = 0.1 GeV, T = 1 GeV")
    plt.legend(fontsize = '20')
    plt.savefig(output_file)

if __name__ == '__main__':
    arg = sys.argv
    plot_spectra_ratio(arg[1], arg[2], arg[3])
