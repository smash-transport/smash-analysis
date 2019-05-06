#/usr/bin/python
# coding=UTF-8
"""
This is a collection of simple functions for analysis
of smash binary output. It is meant to be used by more task-oriented
scripts via:
--------------------------
import read_binary_output
--------------------------
"""

import struct
import sys
import numpy as np
from collections import defaultdict


def _read_binary_header(bfile):
    """Read file header from SMASH binary file."""
    magic, format_version, format_extended, length = struct.unpack('=4sHHi', bfile.read(12))
    if magic != "SMSH":
        print "Fatal error: failed to reproduce magic number."
        sys.exit(1)
    smash_version = struct.unpack('%ds' % length, bfile.read(length))
    assert len(smash_version) == 1
    return smash_version[0], format_extended, format_version


def _read_binary_block_v2(bfile):
    """Read one output block from SMASH binary file."""
    particle_data_type = np.dtype([('p','d',4),('r','d',4),('pdgid','i4'),('id','i4')])

    block_type = bfile.read(1)
    if (block_type == 'p'):
        # got particles block
        npart = np.fromfile(bfile, dtype='i4', count=1)
        particles = np.fromfile(bfile, dtype=particle_data_type, count=npart[0])
        return {'type': block_type,
                'npart': npart[0],
                'part': particles}
    elif (block_type == 'i'):
        # got interaction block
        n_inout = np.fromfile(bfile, dtype='i4', count=2)
        rho = np.fromfile(bfile, dtype='d', count=1)[0]
        sigma = np.fromfile(bfile, dtype='d', count=1)[0]
        incoming = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[0])
        outgoing = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[1])
        return {'type': block_type,
                'nin': n_inout[0],
                'nout': n_inout[1],
                'incoming': incoming,
                'outgoing': outgoing,
                'density': rho,
                'total_cross_section': sigma
               }
    elif (block_type == 'f'):
        n_event = np.fromfile(bfile, dtype='i4', count=1)
        return {'type': block_type,
              'nevent': n_event[0]}
        # got file end block
    elif (block_type == ''):
        # got eof
        return
    else:
        raise ValueError('This is not the start of a block.')

def _read_binary_block_v3(bfile):
    """Read one output block from SMASH binary file."""
    particle_data_type = np.dtype([('p','d',4),('r','d',4),('pdgid','i4'),('id','i4')])

    block_type = bfile.read(1)
    if (block_type == 'p'):
        # got particles block
        npart = np.fromfile(bfile, dtype='i4', count=1)[0]
        particles = np.fromfile(bfile, dtype=particle_data_type, count=npart)
        return {'type': block_type,
                'npart': npart,
                'part': particles}
    elif (block_type == 'i'):
        # got interaction block
        n_inout  = np.fromfile(bfile, dtype='i4', count=2)
        rho      = np.fromfile(bfile, dtype='d',  count=1)[0]
        sigma    = np.fromfile(bfile, dtype='d',  count=1)[0]
        process  = np.fromfile(bfile, dtype='i4', count=1)[0]
        incoming = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[0])
        outgoing = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[1])
        return {'type': block_type,
                'nin': n_inout[0],
                'nout': n_inout[1],
                'incoming': incoming,
                'outgoing': outgoing,
                'density': rho,
                'total_cross_section': sigma,
                'process_type': process}
    elif (block_type == 'f'):
        n_event = np.fromfile(bfile, dtype='i4', count=1)[0]
        return {'type': block_type,
              'nevent': n_event}
        # got file end block
    elif (block_type == ''):
        # got eof
        return
    else:
        raise ValueError('This is not the start of a block.')

def _read_binary_block_v4(bfile):
    """Read one output block from SMASH binary file."""
    particle_data_type = np.dtype([('r','d',4),('mass','d'),('p','d',4),('pdgid','i4'),('id','i4')])

    block_type = bfile.read(1)
    if (block_type == 'p'):
        # got particles block
        npart = np.fromfile(bfile, dtype='i4', count=1)[0]
        particles = np.fromfile(bfile, dtype=particle_data_type, count=npart)
        return {'type': block_type,
                'npart': npart,
                'part': particles}
    elif (block_type == 'i'):
        # got interaction block
        n_inout  = np.fromfile(bfile, dtype='i4', count=2)
        rho      = np.fromfile(bfile, dtype='d',  count=1)[0]
        sigma    = np.fromfile(bfile, dtype='d',  count=1)[0]
        process  = np.fromfile(bfile, dtype='i4', count=1)[0]
        incoming = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[0])
        outgoing = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[1])
        return {'type': block_type,
                'nin': n_inout[0],
                'nout': n_inout[1],
                'incoming': incoming,
                'outgoing': outgoing,
                'density': rho,
                'total_cross_section': sigma,
                'process_type': process}
    elif (block_type == 'f'):
        n_event = np.fromfile(bfile, dtype='i4', count=1)[0]
        return {'type': block_type,
              'nevent': n_event}
        # got file end block
    elif (block_type == ''):
        # got eof
        return
    else:
        raise ValueError('This is not the start of a block.')

def _read_binary_block_v6(bfile):
    """Read one output block from SMASH binary file."""
    particle_data_type = np.dtype([('r','d',4),('mass','d'),('p','d',4),('pdgid','i4'),('id','i4'),('charge','i4')])

    block_type = bfile.read(1)
    if (block_type == 'p'):
        # got particles block
        npart = np.fromfile(bfile, dtype='i4', count=1)[0]
        particles = np.fromfile(bfile, dtype=particle_data_type, count=npart)
        return {'type': block_type,
                'npart': npart,
                'part': particles}
    elif (block_type == 'i'):
        # got interaction block
        n_inout  = np.fromfile(bfile, dtype='i4', count=2)
        rho      = np.fromfile(bfile, dtype='d',  count=1)[0]
        sigma    = np.fromfile(bfile, dtype='d',  count=1)[0]
        sigma_p  = np.fromfile(bfile, dtype='d',  count=1)[0]
        process  = np.fromfile(bfile, dtype='i4', count=1)[0]
        incoming = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[0])
        outgoing = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[1])
        return {'type': block_type,
                'nin': n_inout[0],
                'nout': n_inout[1],
                'incoming': incoming,
                'outgoing': outgoing,
                'density': rho,
                'total_cross_section': sigma,
                'partial_cross_section': sigma_p,
                'process_type': process}
    elif (block_type == 'f'):
        n_event = np.fromfile(bfile, dtype='i4', count=1)[0]
        impact_parameter = np.fromfile(bfile, dtype='d',  count=1)[0]
        return {'type': block_type,
              'nevent': n_event,
              'b' : impact_parameter}
        # got file end block
    elif (block_type == ''):
        # got eof
        return
    else:
        raise ValueError('This is not the start of a block.')

def _read_binary_block_v4_extended(bfile):
    """Read one output block from SMASH binary file."""
    particle_data_type = np.dtype([('r','d',4),('mass','d'),('p','d',4),('pdgid','i4'),('id','i4'),('Ncoll','i4'),('formation_time','d'),('cross_section_scaling_factor','f'),('process_ID_origin','i4'),('process_type_origin', 'i4'),('time_of_origin','f'),('PDG_mother1','i4'),('PDG_mother2','i4')])

    block_type = bfile.read(1)
    if (block_type == 'p'):
        # got particles block
        npart = np.fromfile(bfile, dtype='i4', count=1)[0]
        particles = np.fromfile(bfile, dtype=particle_data_type, count=npart)
        return {'type': block_type,
                'npart': npart,
                'part': particles}
    elif (block_type == 'i'):
        # got interaction block
        n_inout  = np.fromfile(bfile, dtype='i4', count=2)
        rho      = np.fromfile(bfile, dtype='d',  count=1)[0]
        sigma    = np.fromfile(bfile, dtype='d',  count=1)[0]
        process  = np.fromfile(bfile, dtype='i4', count=1)[0]
        incoming = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[0])
        outgoing = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[1])
        return {'type': block_type,
                'nin': n_inout[0],
                'nout': n_inout[1],
                'incoming': incoming,
                'outgoing': outgoing,
                'density': rho,
                'total_cross_section': sigma,
                'process_type': process}
    elif (block_type == 'f'):
        n_event = np.fromfile(bfile, dtype='i4', count=1)[0]
        return {'type': block_type,
              'nevent': n_event}
        # got file end block
    elif (block_type == ''):
        # got eof
        return
    else:
        raise ValueError('This is not the start of a block.')

def _read_binary_block_v5_extended(bfile):
    """Read one output block from SMASH binary file."""
    particle_data_type = np.dtype([('r','d',4),('mass','d'),('p','d',4),('pdgid','i4'),('id','i4'),('Ncoll','i4'),('formation_time','d'),('cross_section_scaling_factor','d'),('process_ID_origin','i4'),('process_type_origin', 'i4'),('time_of_origin','d'),('PDG_mother1','i4'),('PDG_mother2','i4')])

    block_type = bfile.read(1)
    if (block_type == 'p'):
        # got particles block
        npart = np.fromfile(bfile, dtype='i4', count=1)[0]
        particles = np.fromfile(bfile, dtype=particle_data_type, count=npart)
        return {'type': block_type,
                'npart': npart,
                'part': particles}
    elif (block_type == 'i'):
        # got interaction block
        n_inout  = np.fromfile(bfile, dtype='i4', count=2)
        rho      = np.fromfile(bfile, dtype='d',  count=1)[0]
        sigma    = np.fromfile(bfile, dtype='d',  count=1)[0]
        process  = np.fromfile(bfile, dtype='i4', count=1)[0]
        incoming = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[0])
        outgoing = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[1])
        return {'type': block_type,
                'nin': n_inout[0],
                'nout': n_inout[1],
                'incoming': incoming,
                'outgoing': outgoing,
                'density': rho,
                'total_cross_section': sigma,
                'process_type': process}
    elif (block_type == 'f'):
        n_event = np.fromfile(bfile, dtype='i4', count=1)[0]
        return {'type': block_type,
              'nevent': n_event}
        # got file end block
    elif (block_type == ''):
        # got eof
        return
    else:
        raise ValueError('This is not the start of a block.')

def _read_binary_block_v6_extended(bfile):
    """Read one output block from SMASH binary file."""
    particle_data_type = np.dtype([('r','d',4),('mass','d'),('p','d',4),('pdgid','i4'),('id','i4'),('charge','i4'),('Ncoll','i4'),('formation_time','d'),('cross_section_scaling_factor','d'),('process_ID_origin','i4'),('process_type_origin', 'i4'),('time_of_origin','d'),('PDG_mother1','i4'),('PDG_mother2','i4')])

    block_type = bfile.read(1)
    if (block_type == 'p'):
        # got particles block
        npart = np.fromfile(bfile, dtype='i4', count=1)[0]
        particles = np.fromfile(bfile, dtype=particle_data_type, count=npart)
        return {'type': block_type,
                'npart': npart,
                'part': particles}
    elif (block_type == 'i'):
        # got interaction block
        n_inout  = np.fromfile(bfile, dtype='i4', count=2)
        rho      = np.fromfile(bfile, dtype='d',  count=1)[0]
        sigma    = np.fromfile(bfile, dtype='d',  count=1)[0]
        sigma_p  = np.fromfile(bfile, dtype='d',  count=1)[0]
        process  = np.fromfile(bfile, dtype='i4', count=1)[0]
        incoming = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[0])
        outgoing = np.fromfile(bfile, dtype=particle_data_type, count=n_inout[1])
        return {'type': block_type,
                'nin': n_inout[0],
                'nout': n_inout[1],
                'incoming': incoming,
                'outgoing': outgoing,
                'density': rho,
                'total_cross_section': sigma,
                'partial_cross_section': sigma_p,
                'process_type': process}
    elif (block_type == 'f'):
        n_event = np.fromfile(bfile, dtype='i4', count=1)[0]
        impact_parameter = np.fromfile(bfile, dtype='d',  count=1)[0]
        return {'type': block_type,
              'nevent': n_event,
              'b' : impact_parameter}
    elif (block_type == ''):
        # got eof
        return
    else:
        raise ValueError('This is not the start of a block.')

class BinaryReader:
    """A reader for SMASH binary files.

    Use it like this:

    with BinaryReader(path) as reader:
        smash_version = reader.smash_version
        format_version = reader.format_version
        for block in reader:
        ...
    """
    def __init__(self, path):
        self.__file = open(path, 'rb')
        self.smash_version, self.format_extended, self.format_version = _read_binary_header(self.__file)
        if self.format_version == 6:
            if self.format_extended == 1:
                self.__read_block = _read_binary_block_v6_extended
            elif self.format_extended == 0:
                self.__read_block = _read_binary_block_v6
        elif self.format_version == 5:
            if self.format_extended == 1:
                self.__read_block = _read_binary_block_v5_extended
            elif self.format_extended == 0:
                self.__read_block = _read_binary_block_v4
        elif self.format_version == 4:
            if self.format_extended == 1:
                self.__read_block = _read_binary_block_v4_extended
            elif self.format_extended == 0:
                self.__read_block = _read_binary_block_v4
            else:
                print "Fatal error: unknown format variant = ", self.format_extended
                sys.exit(1)
        elif self.format_version == 3:
            self.__read_block = _read_binary_block_v3
        elif self.format_version == 2:
            self.__read_block = _read_binary_block_v2
        else:
            print "Fatal error: unknown format version = ", self.format_version
            sys.exit(1)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.__file.close()

    def __iter__(self):
        return self

    # this is needed for the iterator
    def next(self):
      block = self.read_block()
      if block:
        return block
      else:
        raise StopIteration()

    def read_block(self):
        "Read one output block from SMASH file."
        return self.__read_block(self.__file)


def count_pdg_in_block(block, pdgid):
    """Count particles with given pdgid in a given block."""
    # works only for particles block, not for collisions
    return (block['part']['pdgid'] == pdgid).sum()


def is_elastic22(block):
    """Check if interaction is a 2->2 elastic collision."""
    p_in = block['incoming']['pdgid']
    p_out = block['outgoing']['pdgid']
    return    (  block['nin']  == 2 and \
                 block['nout'] == 2 and \
        ( (p_in[0] == p_out[0] and p_in[1] == p_out[1]) or \
          (p_in[0] == p_out[1] and p_in[1] == p_out[0]) )  )

def is_1to1(block):
    """Check if interaction is 1->1: wall crossing / wall reflection."""
    return (block['nin']  == 1 and block['nout'] == 1)


def get_block_time(block):
    """Get the block time."""
    if (block['type'] == 'p'):
        return block['part']['r'][0][0]
    elif (block['type'] == 'i'):
        return block['incoming']['r'][0][0]
    else:
        print "Error: invalid usage of get_block_time."
        sys.exit(1)

def get_block_E(block):
    """Return energy of particle block."""
    return block['part']['p'][:,0].sum()

def get_block_px(block):
    """Return x component of particle block momentum."""
    return block['part']['p'][:,1].sum()

def get_block_py(block):
    """Return y component of particle block momentum."""
    return block['part']['p'][:,2].sum()

def get_block_pz(block):
    """Return z component of particle block momentum."""
    return block['part']['p'][:,3].sum()

def reaction_Q2(block):
    """Return Q^2 (energy transfer) of interaction."""
    #p0 = block['incoming']['p'][:,0].sum()
    #px = block['incoming']['p'][:,1].sum()
    #py = block['incoming']['p'][:,2].sum()
    #pz = block['incoming']['p'][:,3].sum()
    #return p0*p0 - px*px - py*py - pz*pz

    mom_sqr = np.square(block['incoming']['p'].sum(axis = 0))
    return mom_sqr[0] - mom_sqr[1:].sum()

def reaction_mandelstam_t(block):
    """ Return Mandelstam t of the reaction assuming it is 2->2 reaction. """
    assert(block['nin'] == 2)
    p1 = block['incoming']['p'][0]
    p3 = block['outgoing']['p'][0]
    mom_sqr = np.square(p1 - p3)
    return mom_sqr[0] - mom_sqr[1:].sum()
