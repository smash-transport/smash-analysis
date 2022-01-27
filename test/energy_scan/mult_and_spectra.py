import os
import sys

import numpy as np
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../python_scripts')
import smash_basic_scripts as sb
import argparse
from multiprocessing import Pool


class BulkObservables:
    """
      This class computes and stores some basic bulk observables, obtained from
      SMASH binary output: pt- and rapidity spectra, mean pt at midrapidity,
      multiplicity and midrapidity, mean pt, total multiplicity. All this
      is computed for a given list of hadrons.
    """

    def __init__(self, pdg_list = [211,-211, 111],
                       mtbins = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2]),
                       ybins = np.linspace(-4.0, 4.0, 41),
                       ptbins = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8]),
                       midrapidity_cut = 0.5):
        self.pdglist = pdg_list
        self.pdglist_charha=[211,-211,321,-321,2212,-2212,3122,-3122,1000010020,-1000010020,3312,-3312,3334,-3334,3212,-3212,0]
        self.npdg = len(pdg_list)
        self.mtbins = mtbins
        self.ybins = ybins
        self.ptbins = ptbins
        self.midrapidity_cut = midrapidity_cut
        self.nevents = 0
        self.smash_version = ''

        self.total_multiplicity = np.zeros(self.npdg)
        self.midrapidity_yield = np.zeros(self.npdg)
        self.midrapidity_yield_charha = np.zeros(1)
        self.meanpt_midrapidity = np.zeros(self.npdg)
        self.meanpt_midrapidity_charha = np.zeros(1)
        self.meanmt0_midrapidity = np.zeros(self.npdg)
        self.mthist = np.zeros((self.npdg, mtbins.size - 1))
        self.yhist = np.zeros((self.npdg, ybins.size - 1))
        self.pthist_midrapidity = np.zeros((self.npdg, ptbins.size - 1))
        self.v2 = np.zeros((self.npdg, ptbins.size - 1))

    def __eq__(self, other):
        """ Checks if two sets of bulk observables are identical. Useful for testing. """
        return (self.pdglist == other.pdglist) and \
               abs(self.mtbins - other.mtbins).sum() < 1.e-9 and \
               abs(self.ybins - other.ybins).sum() < 1.e-9 and \
               abs(self.ptbins - other.ptbins).sum() < 1.e-9 and \
               (self.midrapidity_cut == other.midrapidity_cut) and \
               (self.nevents == other.nevents) and \
               (self.smash_version == other.smash_version) and \
               (self.total_multiplicity == other.total_multiplicity).all() and \
               (self.midrapidity_yield == other.midrapidity_yield).all()  and \
               (self.meanpt_midrapidity - other.meanpt_midrapidity < 1.e-5).all() and \
               (self.meanmt0_midrapidity - other.meanmt0_midrapidity < 1.e-5).all() and \
               abs(self.mthist - other.mthist).sum() < 1.e-9 and \
               abs(self.yhist - other.yhist).sum() < 1.e-9 and \
               (self.pthist_midrapidity == other.pthist_midrapidity).all() and \
               (np.abs(self.v2 - other.v2) < 1.e-3).all()

    def updated_mean(self, mean_from_n, mean_from_m, m_over_mplusn):
        """
          Helpful for computing mean of an array A.
          mean_from_n = A[0:n]
          mean_from_m = A[n:n+m]
          m_over_mplusn = m/(n+m)
        """
        return mean_from_n + (mean_from_m - mean_from_n) * m_over_mplusn

    def add_block(self, block):
        """ Read one block of particles and update bulk observables from it. """
        assert(block['type'] == 'p')
        E = block['part']['p'][:,0]
        px = block['part']['p'][:,1]
        py = block['part']['p'][:,2]
        pz = block['part']['p'][:,3]
        m0 = block['part']['mass']
        y = 0.5*np.log((E+pz)/(E-pz))
        pdg = block['part']['pdgid']

        pT2 = px*px + py*py
        pT = np.sqrt(pT2)
        mT0 = np.sqrt(m0*m0 + pT2) - m0
        cos2phi = np.where(pT2 > 0.0, (px*px - py*py) / pT2, 0.0)

        ycut = (np.abs(y) < self.midrapidity_cut)
        self.nevents += 1

        for i in xrange(self.npdg):
            pdgcut = (pdg == self.pdglist[i])
            pdgcut_charha = (pdg == self.pdglist_charha[i])
            added_total_mult = pdgcut.sum()
            added_total_mult_charha = pdgcut_charha.sum()

            if (added_total_mult == 0): continue
            self.yhist[i,:] += np.histogram(y[pdgcut], bins = self.ybins)[0]

            pdg_and_y_cut = np.logical_and(ycut, pdgcut)

            added_midrap_yield = pdg_and_y_cut.sum()


            self.total_multiplicity[i] += added_total_mult
            self.midrapidity_yield[i] += added_midrap_yield


            if (added_midrap_yield == 0): continue
            pt_hist_with_cuts = np.histogram(pT[pdg_and_y_cut], bins = self.ptbins)[0]
            self.mthist[i,:] += np.histogram(mT0[pdg_and_y_cut], bins = self.mtbins)[0]
            self.pthist_midrapidity[i,:] += np.histogram(pT[pdg_and_y_cut], bins = self.ptbins)[0]
            a = float(added_midrap_yield) / self.midrapidity_yield[i]

            self.meanpt_midrapidity[i] = self.updated_mean(self.meanpt_midrapidity[i], pT[pdg_and_y_cut].mean(), a)

            self.meanmt0_midrapidity[i] = self.updated_mean(self.meanmt0_midrapidity[i], mT0[pdg_and_y_cut].mean(), a)

            cos2phi_hist = np.histogram(pT[pdg_and_y_cut], bins = self.ptbins, weights = cos2phi[pdg_and_y_cut])[0]
            self.v2[i,:] += np.where(pt_hist_with_cuts > 0.0, cos2phi_hist, 0.0)

            pdg_and_y_cut_charha = np.logical_and(ycut, pdgcut_charha)
            added_midrap_yield_charha = pdg_and_y_cut_charha.sum()
            if (added_total_mult_charha == 0 or i==self.npdg-1): continue
            self.midrapidity_yield_charha += added_midrap_yield_charha
            b = float(added_midrap_yield_charha) / self.midrapidity_yield_charha
            self.meanpt_midrapidity_charha=self.updated_mean(self.meanpt_midrapidity_charha,pT[pdg_and_y_cut_charha].mean(),b)

    def __iadd__(self, other):
        """ Adds two sets of bulk observables, if it possible. """
        assert(self.pdglist == other.pdglist)
        assert((self.mtbins == other.mtbins).all())
        assert((self.ybins == other.ybins).all())
        assert((self.ptbins == other.ptbins).all())
        assert(self.midrapidity_cut == other.midrapidity_cut)
        assert(self.smash_version == other.smash_version)

        self.nevents += other.nevents
        self.total_multiplicity += other.total_multiplicity
        self.midrapidity_yield += other.midrapidity_yield
        self.mthist += other.mthist
        self.yhist += other.yhist
        self.pthist_midrapidity += other.pthist_midrapidity
        self.v2 += other.v2

        self.meanpt_midrapidity_charha = self.updated_mean(self.meanpt_midrapidity_charha, other.meanpt_midrapidity_charha, a)
        for i in xrange(self.npdg):
            if (other.midrapidity_yield[i] == 0): continue
            a = float(other.midrapidity_yield[i]) / self.midrapidity_yield[i]
            self.meanpt_midrapidity[i] = self.updated_mean(self.meanpt_midrapidity[i], other.meanpt_midrapidity[i], a)
            self.meanmt0_midrapidity[i] = self.updated_mean(self.meanmt0_midrapidity[i], other.meanmt0_midrapidity[i], a)

        return self

    def add_from_file(self, one_file):
        """ Computes bulk observables from a single SMASH output file """
        print one_file
        with sb.BinaryReader(one_file) as reader:
            self.smash_version = reader.smash_version
            block = reader.read_block()
            while block is not None:
                t = block['type']
                assert(t == 'p' or t == 'f')
                if (t == 'p'):
                    # Do not count elastic pp collisions
                    pp_elastic = block['npart'] == 2 and \
                                 (block['part']['pdgid'] == 2212).all()
                    if (not pp_elastic):
                        self.add_block(block)
                block = reader.read_block()

    def add_from_files(self, many_files):
        """ Computes bulk observables from many SMASH output files """
        for f in many_files:
            self.add_from_file(f)

    @staticmethod
    def bin_centers(bin_edges):
        return 0.5 * (bin_edges[1:] + bin_edges[:-1])

    def write_header(self, f):
        f.write('# smash and analysis version: %s %s\n' % \
                (self.smash_version, sb.analysis_version_string()))
        f.write('# total number events: %d\n' % self.nevents)
        f.write('# pdg list:')
        for pdg in self.pdglist: f.write(' %d' % pdg)
        f.write('\n')

    def save(self, files_to_write):
        """ Saves bulk observables to text files. """
        yspectra_file, mtspectra_file, ptspectra_file, v2_file, meanmt0_midrapidity_file,\
        meanpt_midrapidity_file,meanpt_midrapidity_charha_file, midrapidity_yield_file, total_multiplicity_file = files_to_write
        ybin_centers = BulkObservables.bin_centers(self.ybins)
        with open(yspectra_file, 'w') as f:
            self.write_header(f)
            f.write('# y bin edges: %s\n' % np.array_str(self.ybins, max_line_width = 1000000))
            f.write('# y bin center, number of pdg in bin (pdglist)\n')
            for bin_number in xrange(ybin_centers.size):
                f.write(' %6.3f' % ybin_centers[bin_number])
                for i in xrange(self.npdg):
                     f.write(' %8i' % self.yhist[i, bin_number])
                f.write('\n')
        mtbin_centers = BulkObservables.bin_centers(self.mtbins)
        with open(mtspectra_file, 'w') as f:
            self.write_header(f)
            f.write('# mt bin edges: %s\n' % np.array_str(self.mtbins, max_line_width = 1000000))
            f.write('# mt[GeV] bin center, number of pdg in bin (pdglist)\n')
            for bin_number in xrange(mtbin_centers.size):
                f.write(' %6.3f' % mtbin_centers[bin_number])
                for i in xrange(self.npdg):
                    f.write(' %8i' % self.mthist[i, bin_number])
                f.write('\n')
        ptbin_centers = BulkObservables.bin_centers(self.ptbins)
        with open(ptspectra_file, 'w') as f:
            self.write_header(f)
            f.write('# pt bin edges: %s\n' % np.array_str(self.ptbins, max_line_width = 1000000))
            f.write('# pt[GeV] bin center, number of pdg in bin (pdglist)\n')
            for bin_number in xrange(ptbin_centers.size):
                f.write(' %6.3f' % ptbin_centers[bin_number])
                for i in xrange(self.npdg):
                    f.write(' %8i' % self.pthist_midrapidity[i, bin_number])
                f.write('\n')
        with open(v2_file, 'w') as f:
             self.write_header(f)
             f.write('# pt bin edges: %s\n' % np.array_str(self.ptbins, max_line_width = 1000000))
             f.write('# pt[GeV] bin center, v2 of pdg in bin (pdglist)\n')
             for bin_number in xrange(ptbin_centers.size):
                 f.write(' %6.3f' % ptbin_centers[bin_number])
                 for i in xrange(self.npdg):
                     # For each pT, sum of cos2phi over particles has to be divided by number of particles at this pT
                     v2 = np.where(self.pthist_midrapidity[i, bin_number] > 0,
                                   self.v2[i, bin_number] / self.pthist_midrapidity[i, bin_number], 0.0)
                     f.write(' %15.10f' % v2)
                 f.write('\n')

        for i in [meanmt0_midrapidity_file, meanpt_midrapidity_file, meanpt_midrapidity_charha_file,
                  midrapidity_yield_file, total_multiplicity_file]:
            with open(i, 'w') as f: self.write_header(f)
        with open(meanmt0_midrapidity_file, 'a') as f:
            np.savetxt(f, self.meanmt0_midrapidity, fmt = '%8.5f', newline = ' ')
        with open(meanpt_midrapidity_file, 'a') as f:
            np.savetxt(f, self.meanpt_midrapidity, fmt = '%8.5f', newline = ' ')
        with open(meanpt_midrapidity_charha_file, 'a') as f:
            np.savetxt(f, self.meanpt_midrapidity_charha, fmt = '%8.5f', newline = ' ')
        with open(midrapidity_yield_file, 'a') as f:
            np.savetxt(f, self.midrapidity_yield, fmt = '%8i', newline = ' ')
        with open(total_multiplicity_file, 'a') as f:
            np.savetxt(f, self.total_multiplicity, fmt = '%8i', newline = ' ')

    @staticmethod
    def read(files_to_read):
        """ Reads bulk observables, that were saved totext files by the save method. """
        yspectra_file, mtspectra_file, ptspectra_file, v2_file, meanmt0_midrapidity_file,\
        meanpt_midrapidity_file,meanpt_midrapidity_charha_file, midrapidity_yield_file, total_multiplicity_file = files_to_read
        spectra = BulkObservables()
        with open(yspectra_file, 'r') as f:
            smash_version, analysis_version = f.readline().split(':')[1].split()
            spectra.smash_version = smash_version
            spectra.nevents = int(f.readline().split(':')[1])
            spectra.pdglist = [int(x) for x in f.readline().split(':')[1].split()]
            spectra.npdg = len(spectra.pdglist)
            ybinning = f.readline().split('[')[1].split(']')[0].split()
            spectra.ybins = np.array([float(x) for x in ybinning])
            spectra.yhist = np.loadtxt(yspectra_file)[:,1:].T
        with open(mtspectra_file, 'r') as f:
            for _ in xrange(3): f.readline()
            mtbinning = f.readline().split('[')[1].split(']')[0].split()
            spectra.mtbins = np.array([float(x) for x in mtbinning])
            spectra.mthist = np.loadtxt(mtspectra_file)[:,1:].T
        with open(ptspectra_file, 'r') as f:
            for _ in xrange(3): f.readline()
            ptbinning = f.readline().split('[')[1].split(']')[0].split()
            spectra.ptbins = np.array([float(x) for x in ptbinning])
            spectra.pthist_midrapidity = np.loadtxt(ptspectra_file)[:,1:].T
        with open(v2_file, 'r') as f:
            for _ in xrange(3): f.readline()
            f.readline()  # pt-binning is already known from pt-spectra file
            # Multiply back by number of particles at given pT from all events
            spectra.v2 = np.loadtxt(v2_file)[:,1:].T * spectra.pthist_midrapidity

        spectra.meanmt0_midrapidity = np.loadtxt(meanmt0_midrapidity_file)
        spectra.meanpt_midrapidity = np.loadtxt(meanpt_midrapidity_file)
        spectra.meanpt_midrapidity_charha = np.loadtxt(meanpt_midrapidity_charha_file)
        spectra.midrapidity_yield = np.loadtxt(midrapidity_yield_file)
        spectra.total_multiplicity = np.loadtxt(total_multiplicity_file)
        return spectra

    @staticmethod
    def merge_basic_spectra(files_to_read_list, files_to_write):
        """ Merges bulk observables saved in different text files by the save method. """
        spectra = BulkObservables.read(files_to_read_list[0])
        for f in files_to_read_list[1:]:
            spectra += BulkObservables.read(f)
        spectra.save(files_to_write)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--merge", action='store_true',
        help = """
               If not invoked, then input is SMASH binary particle output.
               If invoked, the the input is 8xN files of format 'output_files'
               to be merged into 8 files.
               """)
    parser.add_argument("--output_files", nargs='+', required=True,
        help = """
               Exactly 9 file names in a fixed order:
               output_files = (yspectra_file, mtspectra_file, ptspectra_file, v2_file,
               meanmt0_midrapidity_file,
               meanpt_midrapidity_file,meanpt_midrapidity_charha_file, midrapidity_yield_file, total_multiplicity_file)
               """)
    parser.add_argument("--input_files", nargs='+', required=True,
        help = """
               With --merge these are 8xN filenames of format 'output_files' to be merged.
               Otherwise filename(s) with SMASH binary particles output.
               """)
    parser.add_argument("--parallel", required=False, default=False,
        help = """
               Use paralellized version of the script. Useful if the script is
               run standalone. If run with cmake, then sequential version is enough
               and parallelization is due to cmake
               """)
    args = parser.parse_args()

    if (args.merge):
        assert(len(args.input_files) % 8 == 0)
        files_to_read_list = zip(*(iter(args.input_files),) * 8)
        BulkObservables.merge_basic_spectra(files_to_read_list, args.output_files)
    else:
        pdg_list = [211,-211,111,321,-321,2212,-2212,3122,-3122,1000010020,-1000010020,3312,-3312,3334,-3334,3212,-3212]

        def get_bulk_observables_from_file(input_file):
            b = BulkObservables(pdg_list = pdg_list)
            b.add_from_file(input_file)
            return b
        if (args.parallel):
           pool = Pool()
           spectra = get_bulk_observables_from_file(args.input_files[0])
           spectra_list = pool.map_async(get_bulk_observables_from_file, args.input_files[1:])
           for i in spectra_list.get():
               spectra += i
        else:
            spectra = BulkObservables(pdg_list = pdg_list)
            spectra.add_from_files(args.input_files)
        spectra.save(args.output_files)

        # Read-write consistency check
        spectra2 = BulkObservables.read(args.output_files)
        assert(spectra == spectra2)
