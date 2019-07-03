#!/usr/bin/python2
# -*- encoding: utf-8 -*-

"""
geometrical_cross_section.py

This script reads SMASH binary output
and prints out the cross sections for processes
which have the given PDG code(s) in initial state.

SMASH uses the geometrical interpretation of cross section
to calculate the interaction distance between two particles.
Here the procedure is "reverse-engineered"
to extract the cross section from the output data:
1. Boost particles to center-of-momentum frame
2. Calculate transverse distance R
3. Add the transverse distance R to the event histogram
4. Add the transverse distance R to the collision histogram if
   a collision takes place
5. Get the probability of scattering F(R) = scattering list / event list
6. Calculate the expectation value of the total cross section as
   <S> = - pi \int R^2 dF(R),
   and the square of its fluctuation as
   <\Delta S^2> = <S^2> - <S>^2, where
   <S^2> = - pi^2 \int R^4 dF(R)
   (See issue #5799)


Usage: ./geometrical_cross_section.py [options] [pdg code(s)] [output file(s)]
Example: p+p: ./geometrical_cross_section.py 2212 2212 data/0/collisions_binary.bin

Options:
--verbose       Print additional information in special situations

Notice:
This script only calculates the cross section with a certain colliding
energy. If you'd like to scan the energies within a certain region, please
use "make" to build the cross section target, where the cross sections are
calculated for different colliding energies and finally combined using a
combining script.
"""

import sys
import os
import argparse
import math
import numpy
from collections import defaultdict, Counter, deque
sys.path.append(os.path.dirname(
    os.path.abspath(__file__)) + '/../../python_scripts')
import smash_basic_scripts as smash
import reconstruct_interaction_graph as ig
from lorentz import lorentz_boost
from pdgs_from_config import charge_str_neg, charge_str_zero, charge_str_pos

reload(sys)
sys.setdefaultencoding('utf-8')

numpy.seterr(all='raise')

bmax = 2.5
nbins = 1000
b_bins = numpy.linspace(0, bmax, nbins+1)[1:]
b_bins_tot = numpy.zeros(nbins)
b_bins_scat = numpy.zeros(nbins)

def transverse_distance(
        momentum_a, momentum_b, position_a, position_b):
    """Compute transverse distance between particles in their local rest frame.

    This function computes the transverse distance between two particles
    as defined in UrQMD (See arXiv:nucl-th/9803035, Eq. (3.27)).

    \param[in] momentum_a List containing 4-momentum components of particle a.
    \param[in] momentum_b List containing 4-momentum components of particle b.
    \param[in] position_a List containing spacetime position of particle a.
    \param[in] position_b List containing spacetime position of particle b.

    \return Transverse distance between a and b in center-of-momentum frame.
    """

    # compute center-of-momentum velocity
    velocity_cm = ((momentum_a[1:4] + momentum_b[1:4]) /
                   (momentum_a[0] + momentum_b[0]))

    # Lorentz boost to center-of-momentum frame (if not there already)
    if abs(velocity_cm).any() > 0.0001:
        momentum_a = lorentz_boost(momentum_a, velocity_cm)
        momentum_b = lorentz_boost(momentum_b, velocity_cm)
        position_a = lorentz_boost(position_a, velocity_cm)
        position_b = lorentz_boost(position_b, velocity_cm)

    # three-vectors in CM frame
    p_a = numpy.array(momentum_a[1:4])
    p_b = numpy.array(momentum_b[1:4])
    x_a = numpy.array(position_a[1:4])
    x_b = numpy.array(position_b[1:4])

    try:
        result = (numpy.inner(x_a - x_b, x_a - x_b) -
                  (numpy.inner(x_a - x_b, p_a - p_b))**2 /
                  numpy.inner(p_a - p_b, p_a - p_b))
    except FloatingPointError:
        print >> out, "FloatingPointError! p_a:", p_a, "p_b:", p_b
        print >> out, "velocity:", velocity_cm
        result = 0
    result = max(0, result)

    return math.sqrt(result)


def parse_arguments():
    """Parse and return the command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="print additional output")
    parser.add_argument("-f", "--first", action="store_true",
                        help="count only the first interaction in an event")
    parser.add_argument(
        "-o1",
        "--output1",
        type=str,
        help="set first output file: individual channels (instead of printing to stdout)")
    parser.add_argument(
        "-o2",
        "--output2",
        type=str,
        help="set second output file: grouped channels (instead of printing to stdout)")
    parser.add_argument(
        "-o3",
        "--output3",
        type=str,
        help="set third output file: individual final channels (instead of printing to stdout)")
    parser.add_argument(
        "-o4",
        "--output4",
        type=str,
        help="set fourth output file: grouped final channels (instead of printing to stdout)")
    parser.add_argument(
        "-o5",
        "--output5",
        type=str,
        help="set fifth output file: individual final channels for unstable final particles (instead of printing to stdout)")
    parser.add_argument(
        "-o6",
        "--output6",
        type=str,
        help="set sixth output file: grouped final channels for unstable final particles (instead of printing to stdout)")
    parser.add_argument(
        "-o7",
        "--output7",
        type=str,
        help="set seventh output file: process types (instead of printing to stdout)")
    parser.add_argument("-r", "--resonances", help="count reconstructed resonances", type=str)
    parser.add_argument("pdg1", help="PDG code of the first particle")
    parser.add_argument("pdg2", help="PDG code of the second particle")
    parser.add_argument("config_file", help="config file")
    parser.add_argument("filename", nargs='+',
                        help="binary file(s) containing collision history")
    args = parser.parse_args()
    return args


def pdgs_to_names(pdgs, pdg_to_name):
    """Translate the given PDG codes to particles names."""
    names = sorted([pdg_to_name(pdg) for pdg in pdgs])
    return names


def get_energyindex(energy):
    """Calculate the energy index given an energy in GeV.

    This is important, since the numbers of reactions are counted per energy.
    Small difference between the final and initial state energies due to
    numerics could result the multiplicity not summing up correctly, leading to
    artificial dips in the final state cross sections.
    """
    # Round based on SMASH's tolerance for energy conservation violation
    return round(energy, 5)


def get_sqrts(ptot):
    """Calculate the center-of-mass energy given the total momentum."""
    return math.sqrt(ptot[0]**2 - ptot[1]**2 - ptot[2]**2 - ptot[3]**2)


def get_distance_index(datablock):
    """Calculate the bin index of the transverse distance."""
    # extract the 4-vectors of initial particles
    if datablock['type'] == 'i':
       momentum_a = datablock['incoming']['p'][0, :]
       momentum_b = datablock['incoming']['p'][1, :]
       position_a = datablock['incoming']['r'][0, :]
       position_b = datablock['incoming']['r'][1, :]
    elif datablock['type'] == 'p':
       momentum_a = datablock['part']['p'][0, :]
       momentum_b = datablock['part']['p'][1, :]
       position_a = datablock['part']['r'][0, :]
       position_b = datablock['part']['r'][1, :]

    ptot = momentum_a + momentum_b

    sqrts = get_sqrts(ptot)
    energyindex = get_energyindex(sqrts)

    # compute distance between initial particles
    d_T = transverse_distance(momentum_a, momentum_b,
                              position_a, position_b)
    id_T = int(d_T // (bmax / nbins))

    return id_T, energyindex


def is_first_reaction(datablock):
    """Does the given data block correspond to the first reaction?"""
    return (datablock['incoming']['id'] == [0, 1]).all()


def update_process_list(process_list, tag, pdg_to_name, final_pdgs):
    """Add process name to list and return it."""
    final_particles = pdgs_to_names(final_pdgs, pdg_to_name)
    process = "+".join(final_particles)
    process_list[tag].add(process)
    return process

def dict_to_str(d):
    """Convert a dict to a human-readable string."""
    return '{' + ',\n'.join('{}: {}'.format(k, v) for k, v in sorted(d.iteritems())) + '}'

def list_to_str(l):
    """Convert a list to a human-readable string."""
    return '[' + ', '.join(str(i) for i in l) + ']'

def counter_to_str(c):
    """Convert a counter to a string by repeating."""
    return '+'.join('+'.join([key] * count) for key, count in sorted(c.iteritems()))

def calc_tot_xs():
    """
     Calculate the expectation value and fluctuation oftotal cross section using
     <S> = - pi \int R^2 dF(R),
     and the square of its fluctuation as
     <\Delta S^2> = <S^2> - <S>^2, where
     <S^2> = - pi^2 \int R^4 dF(R)
     (See issue #5799)
    """
    b_bins_prob = b_bins_scat[b_bins_tot != 0] / b_bins_tot[b_bins_tot != 0]
    xs_tot = - numpy.sum(b_bins[b_bins_tot != 0][:-1] ** 2 * numpy.diff(b_bins_prob)) * math.pi * 10
    xs_tot_2 = - numpy.sum(b_bins[b_bins_tot != 0][:-1] ** 4 * numpy.diff(b_bins_prob)) * (math.pi * 10) ** 2
    d_xs_tot_2 = xs_tot_2 - xs_tot ** 2
    # Due to numeric inaccuracy, the fluctuation square can sometimes be a very small negative number
    # (such as 10^{-14}), so it's set equal to zero in this case.
    d_xs_tot_2 = max(0.0, d_xs_tot_2)
    return xs_tot, d_xs_tot_2

class Processor:
    """Store the variables required for processing events to calculate cross
    sections.

    This has the advantage the all the methods can use these variables without
    having to pass them as arguments.
    """
    def __init__(self, resonances=[]):
        self.args = parse_arguments()
        self.process_list = defaultdict(set)  # list of processes for each tag
        # Currently used tags are:
        # - 'individual': channels, with resonances
        # - 'grouped': generic channels, with resonances
        # - 'final_individual': final-state channels, only stable particles
        # - 'final_grouped': generic final-state channels, only stable particles
        # - 'total': total count of processes with the initial pdgs
        # - 'elastic': number of elastic collisions
        # number of times a process of given tag has happened
        # Currently used tags are:
        # - 'elastic'
        # - 'via resonance'
        # - 'parametrized inelastic'
        # - 'soft string'
        # - 'hard string'
        # types such as 'nothing', 'resonance decay', 'wall transition' are not
        # included since they don't happen initially.
        self.energyindex = 0
        # This is used to check whether all the events have the same energy
        self.energyindex_old = None
        self.process_count = defaultdict(int)
        self.type_count = defaultdict(int)
        # list of pdg codes, read as input
        self.initial_pdgs = sorted([int(self.args.pdg1), int(self.args.pdg2)])
        # list of particle names
        self.initial_parts = sorted(map(self.get_name, self.initial_pdgs))
        self.done = False

        # graph representing all particles and their interactions
        # (filled in `count_interaction`)
        self.interaction_graph = defaultdict(ig.Particle)
        resonances = self.args.resonances.split(',')

        def str_to_unicode(s):
            return unicode(s, 'utf-8')

        self.r_unicodes = [str_to_unicode(r) for r in resonances]

        self.observable_resonances = 0
        self.unobservable_resonances = 0

    def get_name(self, pdg):
        """Get the name of a particle given the PDG ID."""
        return smash.pdg_to_name(pdg, self.args.config_file)

    def get_name_generic(self, pdg):
        """Get the generic name of a particle given the PDG ID."""
        return smash.pdg_to_name_generic(pdg, self.args.config_file)

    def generic_name(self, pdgid):
        """Get the generic name with charge.

        For example, instead of 'N(1900)', write 'N*'
        """
        return smash.generic_name(self.get_name(pdgid))

    def unify_name(self, pdgid):
        """Get the generic name without charge."""
        return smash.strip_charge(self.generic_name(pdgid))

    def count_initial_interaction(self, datablock):
        """Count the processes in the initial interaction."""
        id_T, self.energyindex = get_distance_index(datablock)
        # Check whether the event has the same energy as the previous event
        assert (self.energyindex == self.energyindex_old or self.energyindex_old == None)
        self.energyindex_old = self.energyindex

        if datablock['type'] == 'i':
           b_bins_scat[id_T] += 1
           final_pdgs = datablock['outgoing']['pdgid']
           process = update_process_list(
               self.process_list, 'individual', self.get_name, final_pdgs)
           process_gen = update_process_list(
               self.process_list, 'grouped', self.unify_name, final_pdgs)
           self.process_name = process
           self.process_name_generic = process_gen
           self.process_count[('individual', process)] += 1
           self.process_count[('grouped', process_gen)] += 1
           self.process_count['total'] += 1
           type_id = datablock['process_type']
           # Get the type of the action. Elastic scatterings are not counted
           # here, they'll be counted according to the final states
           if (type_id == 1):
              p_type = 'elastic'
           elif (type_id == 2 or type_id == 3):
              # '2' corresponds to a 2->1 process, '3' to a 2->2 process, both inelastic
              # Now we make the following assumption:
              # We expect, that the in both cases, at least one resonance is
              # formed, so that we denote it 'resonance' in the legend.
              # This is justified from the low-energy contributions to the
              # cross sections coming from resonance formations and decays (in
              # pp, pn, pipi, pip, ... processes) only. This assumption
              # will however not hold anymore if there will be processes like
              # 'p + n -> p + n + pi' or similar. But as there are none implemented
              # in SMASH at the moment, we assume the above assumption is justified.
              #
              # In the case of the kaon-nucleon interactions, the story is different:
              # Since the inelastic 2->2 K-N cross sections are parameterized,
              # they do not rely purely on resonance formations.
              # Here, we will denote those processes as 'parametrized-inelastic' to
              # make it clear, we are using a parametrization.
              # In addition we need to trigger on the inital state particles to
              # see if there are kaons involved and whether we are dealing
              # with  inelastic 2->2 processes that are unparametrized or
              # parametrized.

              initial_state_PDGs = [datablock['incoming'][0][3], datablock['incoming'][1][3]]
              if (321 in initial_state_PDGs or -321 in initial_state_PDGs): # 321 corresponds to kaons
                  p_type = 'parametrized-inelastic'
              else:
                  p_type = 'resonance'
           elif (type_id >= 41 and type_id <= 45):
              p_type = 'soft-string'
           elif (type_id == 46):
              p_type = 'hard-string'
           if p_type != 'elastic':
              self.type_count[p_type] += 1

           # determine masses of initial particles
           self.initial_masses = map(lambda x: x[1],
                                     sorted((p['id'], p['mass']) for p in datablock['incoming']))
        elif datablock['type'] == 'p':
           b_bins_tot[id_T] += 1
           self.initial_masses = map(lambda x: x[1],
                                     sorted((p['id'], p['mass']) for p in datablock['part']))

    def count_final_state(self, datablock):
        """Count how often this final state occurs.

        Updates the corresponding process count for 'elastic' and 'final_grouped'."""
        assert datablock['type'] == 'p'
        assert isinstance(self.initial_pdgs, list)
        final_pdgs = sorted(datablock['part']['pdgid'])

        # determine total 4-mom
        ptot = datablock['part']['p'].sum(axis=0)
        sqrts = get_sqrts(ptot)

        process_individual = update_process_list(
            self.process_list,
            'final_individual',
            self.get_name,
            final_pdgs)
        process_generic = update_process_list(
            self.process_list,
            'final_grouped',
            self.unify_name,
            final_pdgs)

        # elastic and total cross section
        if datablock['npart'] == 2 and self.initial_pdgs == final_pdgs:
            # momentum of first particle
            momentum_a_out = datablock['part']['p'][0, :]
            # only count particles with transverse momentum
            if abs(momentum_a_out[1]) > 1e-6 or abs(momentum_a_out[2]) > 1e-6:
                self.process_count['elastic'] += 1
                self.type_count['elastic'] += 1
                self.process_count[('final_individual',
                                    process_individual)] += 1
                self.process_count[(
                    'final_grouped', process_generic)] += 1
        else:
            self.process_count[('final_individual',
                                process_individual)] += 1
            process_individual_per_res = '{}→{}'.format(self.process_name, process_individual)
            self.process_list['final_individual'].add(process_individual_per_res)
            self.process_count[('final_individual', process_individual_per_res)] += 1

            self.process_count[(
                'final_grouped', process_generic)] += 1
            process_generic_per_res = '{}→{}'.format(smash.strip_charge(self.process_name_generic), process_generic)
            self.process_list['final_grouped'].add(process_generic_per_res)
            self.process_count[('final_grouped', process_generic_per_res)] += 1

        # unstable final states

        final_individual = Counter(map(self.get_name, final_pdgs))
        final_generic = Counter(map(self.unify_name, final_pdgs))

        # Calculate the maximal number of mothers in the reaction chain from
        # final to initial state.
        max_indeg = ig.max_indegrees_bottom_up(
            datablock['part']['id'], self.interaction_graph)

        for pid, particle in self.interaction_graph.iteritems():
            go_on = True
            for r in self.r_unicodes:
                if go_on and smash.strip_charge(self.get_name_generic(particle.pdgid)) == r:
                    go_on = False
                    is_observable = all(max_indeg[d] <= 1 for d in set(particle.daughters) - set(particle.mothers))
                    if is_observable:
                        self.observable_resonances += 1
                    else:
                        self.unobservable_resonances += 1
                        # We don't count unobservable resonances.
                        # Also, we can't substitute the daughters by the resonance
                        # if there are 2 -> 2 reactions involved!
                        continue
                    '''
                    print >> sys.stderr, dict_to_str(self.interaction_graph)
                    print >> sys.stderr, dict_to_str(max_indeg)
                    print >> sys.stderr, 'is observable:', is_observable
                    '''

                    daughters = ig.get_daughters(self.interaction_graph, pid)
                    daughters_count = Counter(self.get_name(self.interaction_graph[p].pdgid) for p in daughters)
                    motherid = pid


                    individual_mother_name = self.get_name(self.interaction_graph[motherid].pdgid)
                    individual_daughter_names = [self.get_name(self.interaction_graph[i].pdgid) for i in daughters]
                    process_res_individual = counter_to_str(final_individual - Counter(individual_daughter_names) + Counter([individual_mother_name]))
                    assert r in process_res_individual, '"{}" not in "{}"'.format(r, process_res_individual)
                    # It seems that only the symbol of charge is removed from the process_res_individual string. Not sure
                    # whether it's good to call it "generic".
                    process_res_generic = (process_res_individual
                        .replace(charge_str_pos, '')
                        .replace(charge_str_zero, '')
                        .replace(charge_str_neg, '')
                    )
                    assert r in process_res_generic, '"{}" not in "{}"'.format(r, process_res_generic)

                    self.process_list['final_individual_res'].add(process_res_individual)
                    self.process_count[('final_individual_res', process_res_individual)] += 1
                    process_res_individual_per_res = '{}→{}'.format(self.process_name, process_res_individual)
                    self.process_list['final_individual_res'].add(process_res_individual_per_res)
                    self.process_count[('final_individual_res', process_res_individual_per_res)] += 1

                    self.process_list['final_grouped_res'].add(process_res_generic)
                    self.process_count[('final_grouped_res', process_res_generic)] += 1

                    process_res_generic_per_res = '{}→{}'.format(smash.strip_charge(self.process_name_generic), process_res_generic)
                    self.process_list['final_grouped_res'].add(process_res_generic_per_res)
                    self.process_count[('final_grouped_res', process_res_generic_per_res)] += 1
                    # We want to consider the special case where a phi was
                    # reconstructed and we are interested in the final state.
                    # This is important for comparison to experiment, because
                    # phis are usually only reconstructed from K+ K- only.
                    if process_res_individual == u'N⁺+N⁺+φ':
                        process_phi_individual = u'{}→φ+N⁺+N⁺→{}'.format(
                                self.process_name, process_individual)
                        process_phi_grouped = u'{}→φ+N⁺+N⁺→{}'.format(
                                smash.strip_charge(self.process_name_generic),
                                process_generic)
                        self.process_list['final_individual_res'].add(process_phi_individual)
                        self.process_count[('final_individual_res', process_phi_individual)] += 1
                        self.process_list['final_grouped_res'].add(process_phi_grouped)
                        self.process_count[('final_grouped_res', process_phi_grouped)] += 1

    def process_datablocks(self, datablocks):
        """Process all data blocks and update `process_count`."""
        reading_allowed = True
        seen_particle_block = 0
        event = 0
        seen_interaction = False
        empty_events = 0

        for i, datablock in enumerate(datablocks):
            if datablock['type'] == 'i' and reading_allowed:
                if is_first_reaction(datablock):
                    self.count_initial_interaction(datablock)
                    if self.args.first:
                        reading_allowed = False
                ig.count_interaction(self.interaction_graph, datablock)
                seen_interaction = True
            elif datablock['type'] == 'p':
                assert seen_particle_block < 2
                # the first block are the initial particles,
                # the second block the final particles
                if seen_particle_block == 0:
                    # bin the transverse distance
                    # determine the types (strings, resonances, etc.) of scattering
                    self.count_initial_interaction(datablock)
                elif seen_particle_block == 1:
                    # determine elastic cross section from final particle list,
                    # where all resonances have been forced to decay
                    self.count_final_state(datablock)
                seen_particle_block += 1
            elif datablock['type'] == 'f':
                reading_allowed = True
                seen_particle_block = 0
                event += 1
                if not seen_interaction:
                    empty_events += 1
                seen_interaction = False
                self.interaction_graph = defaultdict(ig.Particle)
        print '{:.1f}% = {} / {} of the events were empty.'.format(
                float(empty_events) / event * 100, empty_events, event)

    def print_to_file(self, out, tag, process_el=""):
        """Print data (process_list, process_count) to output file 'out'."""
        process_list = sorted(self.process_list[tag])
        # look for the possible elastic process and move it to the front
        if process_el:
            try:
                index_elast = process_list.index(process_el)
                process_list.pop(index_elast)
                process_list.insert(0, process_el)
            except IndexError:
                if self.args.verbose:
                    print >> sout, "# Note: No elastic process found."
            except ValueError:
                if self.args.verbose:
                    print >> sout, "# Note:", process_el, "(elastic process) is not in process list."
                    print >> sout, "# Following processes were found:",
                    for item in process_list:
                        print >> sout, item,
                    print >> sout
        # write initial state
        print >> out, "#initial",
        for pdg in self.initial_parts:
            print >> out, pdg,
        print >> out
        # write initial masses
        print >> out, "#masses",
        assert len(self.initial_masses) == 2
        print >> out, ' '.join(str(mass) for mass in self.initial_masses)
        # write version
        print >> out, "#version", self.smashversion, "format", self.formatversion, \
            "analysis", smash.analysis_version_string()
        # write labels
        print >> out, "#columns $\\sqrt{s}$[GeV] total total_err",
        if tag == 'process_type':
           for p_type in sorted(self.type_count.iterkeys()):
               print >> out, p_type, p_type + '_err',
        else:
           if process_el and (process_el not in process_list):
               print >> out, process_el, process_el + '_err',   # elastic
           for proc in process_list:
               print >> out, proc, proc + '_err',
        print >> out
        # write cross sections and errors
        total_count = self.process_count['total']
        if total_count > 0:
          total_xsection, total_xsection_err_sq = calc_tot_xs()
          #assert total_xsection_err_sq >= 0.
          total_xsection_rel_err_sq = total_xsection_err_sq / (total_count ** 2)
          total_xsection_err = math.sqrt(total_xsection_err_sq)
          print >> out, self.energyindex, total_xsection, total_xsection_err,
          if tag == 'process_type':
             for p_type in sorted(self.type_count.iterkeys()):
                count = float(self.type_count[p_type])
                xs = total_xsection * count / total_count
                excl_xsection_rel_err_sq = 1./count if count else 0.
                error = xs * numpy.sqrt( total_xsection_rel_err_sq + excl_xsection_rel_err_sq )
                print >> out, xs, error,
          else:
             if process_el and (process_el not in process_list):
                count = float(self.process_count['elastic'])
                xs = total_xsection * count / total_count
                excl_xsection_rel_err_sq = 1./count if count else 0.
                error = xs * numpy.sqrt( total_xsection_rel_err_sq + excl_xsection_rel_err_sq )
                print >> out, xs, error,
             for proc in process_list:
                count = float(self.process_count[(tag, proc)])
                xs = total_xsection * count / total_count
                excl_xsection_rel_err_sq = 1./count if count else 0.
                error = xs * numpy.sqrt( total_xsection_rel_err_sq + excl_xsection_rel_err_sq )
                print >> out, xs, error,

          if self.args.verbose:
              if total_xsection < 0.1:
                  print >> sys.stderr, "Warning: low cross section encountered: energy: {}, cross section: {}".format(
                      self.energyindex, total_xsection)
          print >> out

    def print_unobservable_resonances(self):
        """Print how many resonances were not observable."""
        total_resonances = self.observable_resonances + self.unobservable_resonances
        if total_resonances:
            print 'unobservable resonances: {}/{} = {:.1f}%'.format(
                self.unobservable_resonances, self.observable_resonances,
                float(self.unobservable_resonances) / total_resonances * 100)

    def process_all(self):
        """Process all input files and output counts."""
        if self.done:
            raise RuntimeError(
                '`Process.process_all` should be called only once')
        firstfile = True
        reading_allowed = True
        # go through the list of input arguments
        for name in self.args.filename:
            # if the argument is a file, read the contents
            if os.path.isfile(name):
                with smash.BinaryReader(name) as reader:
                    # check the header for version information
                    self.smashversion = reader.smash_version
                    self.formatversion = reader.format_version
                    if firstfile:
                        smash_previous = self.smashversion
                        format_previous = self.formatversion
                        firstfile = False
                    elif (self.smashversion != smash_previous or
                          self.formatversion != format_previous):
                        if self.args.verbose:
                            print >> sout, "# Data from different versions detected!"
                            print >> sout, "# Above version from file", name
                    # loop over all data blocks in the file
                    self.process_datablocks(reader)
                self.print_unobservable_resonances()
        process_elast = "+".join(self.initial_parts)
        process_elast_gen = "+".join(pdgs_to_names(self.initial_pdgs,
                                                   self.unify_name))

        # set up output streams
        out1 = open(
            self.args.output1,
            'w') if self.args.output1 else sys.stdout
        out2 = open(
            self.args.output2,
            'w') if self.args.output2 else sys.stdout
        out3 = open(
            self.args.output3,
            'w') if self.args.output3 else sys.stdout
        out4 = open(
            self.args.output4,
            'w') if self.args.output4 else sys.stdout
        out5 = open(
            self.args.output5,
            'w') if self.args.output5 else sys.stdout
        out6 = open(
            self.args.output6,
            'w') if self.args.output4 else sys.stdout
        out7 = open(
            self.args.output7,
            'w') if self.args.output5 else sys.stdout

        self.print_to_file(out1, 'individual', 'elastic')
        self.print_to_file(out2, 'grouped', 'elastic')
        self.print_to_file(out3, 'final_individual', 'elastic')
        self.print_to_file(out4, 'final_grouped', 'elastic')
        self.print_to_file(out5, 'final_individual_res', 'elastic')
        self.print_to_file(out6, 'final_grouped_res', 'elastic')
        self.print_to_file(out7, 'process_type')

        self.done = True


if __name__ == '__main__':

    sout = sys.stdout

    p = Processor()
    p.process_all()
