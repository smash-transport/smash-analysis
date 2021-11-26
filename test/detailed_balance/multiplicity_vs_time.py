"""
Purpose: get multiplicity of every particle in a given list as a function of
         time from multiple SMASH output files.
"""

import sys
import os  # operating system interface
import numpy as np
import argparse  # command line argument parser
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts')
import smash_basic_scripts as sb
from multiprocessing import Pool

parser = argparse.ArgumentParser()
# options and arguments
parser.add_argument("output_file", type=str)
parser.add_argument("pdg_list", type=str)
parser.add_argument("config_file", help="config file")
parser.add_argument("files_to_analyze", nargs='+',
                    help="binary file(s) containing collision history")
args = parser.parse_args()

pdg_list = np.array([int(sb.name_to_pdg(x, args.config_file)) for x in args.pdg_list.split(',')])
total_pdgs = pdg_list.shape[0]
max_tcounter = 10000  # max number of blocks per event

def get_multiplicity(file_to_analyze):
    """
        Compute number of particles against time in the file file_to_analyze
        for every particle in a pdg_list.
    """
    time    = np.zeros(max_tcounter)
    mul     = np.zeros((max_tcounter, total_pdgs))
    mul_sqr = np.zeros((max_tcounter, total_pdgs))
    event_num = 0
    tcounter = 0
    with sb.BinaryReader(file_to_analyze) as reader:
        smash_version = reader.smash_version
        #format_version = reader.format_version
        for block in reader:
            if (block['type'] == 'f'):  # end of event
                event_num += 1
                blocks_per_event = tcounter
                tcounter = 0
            if (block['type'] == 'i'):  # interaction
                print('Error: there should be no interactions in this file!')
                sys.exit(1)
            if (block['type'] == 'p'):  # particles
                if (event_num == 0):
                    time[tcounter] = sb.get_block_time(block)
                for i in np.arange(total_pdgs):
                    part_num = sb.count_pdg_in_block(block, pdg_list[i])
                    mul[tcounter, i] += part_num
                    mul_sqr[tcounter, i] += part_num*part_num
                tcounter += 1
    return { 'smash_version' : smash_version,
             'nevents' : event_num,
             'blocks_per_event' : blocks_per_event,
             'time' : time,
             'mul' : mul,
             'mul_sqr' : mul_sqr
           }

# Non-paralellized alternatives:
# 1.
#  results = map(get_multiplicity, args.files_to_analyze)
# 2.
#  results = []
#  for file_to_analyze in args.files_to_analyze:
#      results.append(get_multiplicity(file_to_analyze))
pool = Pool()
results = pool.map_async(get_multiplicity, args.files_to_analyze)

event_num = 0
blocks_per_event = 0
time    = np.zeros(max_tcounter)
mul     = np.zeros((max_tcounter, total_pdgs))
mul_sqr = np.zeros((max_tcounter, total_pdgs))

for res in results.get():
    event_num += res['nevents']
    if (blocks_per_event == 0):
        blocks_per_event = res['blocks_per_event']
        time = res['time']
        smash_version = res['smash_version']
    else:
        assert(blocks_per_event == res['blocks_per_event'])
        assert(smash_version == res['smash_version'])
    mul += res['mul']
    mul_sqr += res['mul_sqr']

mul /= event_num  # get average multiplicity
mul_sqr /= event_num # get average square of multiplicity
mul_err = np.sqrt((mul_sqr - mul * mul)/(event_num - 1))
with open(args.output_file, 'w') as f:
    f.write('# smash and analysis version\n')
    f.write('%s %s\n' % (smash_version, sb.analysis_version_string()))
    f.write('# total number events\n')
    f.write('%d\n' % event_num)
    f.write('# pdg codes list\n')
    for i in np.arange(total_pdgs): f.write('%d ' % pdg_list[i])
    f.write('\n')
    f.write('# time moments array [%d]\n' % blocks_per_event)
    for i in np.arange(blocks_per_event): f.write('%.3f ' % time[i])
    f.write('\n')
    for i in np.arange(total_pdgs):
        f.write('# multiplicities and their stat errors of pdg %d versus time\n' % pdg_list[i])
        for j in np.arange(blocks_per_event): f.write('%.5f ' % mul[j, i])
        f.write('\n')
        for j in np.arange(blocks_per_event): f.write('%.5f ' % mul_err[j, i])
        f.write('\n')
