import os
import sys
import yaml
import fnmatch
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts')
import smash_basic_scripts as sb
import numpy as np

def count_elastic_scat(file_to_analyze, t_start):
    with sb.BinaryReader(file_to_analyze) as reader:
        smash_version = reader.smash_version
        format_version = reader.format_version
        elastic = 0
        for block in reader:
            if (block['type'] == 'f'):  # end of event
                event_num = block['nevent'] + 1
                # print "finished analyzing event ", event_num
            if (block['type'] == 'i' and sb.is_elastic22(block)): # 2->2 el. scatt.
                time = sb.get_block_time(block)
                if (time > t_start): elastic += 1
    return { 'smash_ver' : smash_version,
             'ev_num' : event_num,
             'elast_coll' : elastic}

# Get arguments: parent directory of smash output, output file name
data_directory = sys.argv[1]
output_file    = sys.argv[2]

# To avoid 'warm-up' effects, start counting after some time passes
t_start = 5.0  # fm/c

# Find all collisions_binary.bin files in the folder
paths_to_analyze = []
for root, dirnames, filenames in os.walk(data_directory):
  for filename in fnmatch.filter(filenames, 'collisions_binary.bin'):
      paths_to_analyze.append(root)

outfile = open(output_file,'w')
outfile.write('# Ncoll(2->2 elastic) N_ev run_time[fm/c] V[fm^3]'
              ' sigma[fm^2] N_part  N_test Temperature[GeV] timestep[fm/c]\n')

assert paths_to_analyze
for path in paths_to_analyze:
    # Get the dependency variable from config
    conf_file = os.path.join(path,'config.yaml')
    with open(conf_file, 'r') as cf:
        smash_conf = yaml.load(cf, Loader=yaml.loader.BaseLoader)
    L = float(smash_conf["Modi"]["Box"]["Length"])
    V = L*L*L
    N_part = int(smash_conf["Modi"]["Box"]["Init_Multiplicities"]["111"])
    sigma  = float(smash_conf["Collision_Term"]["Elastic_Cross_Section"]) * 0.1  # mb -> fm^2
    N_test = int(smash_conf["General"]["Testparticles"])
    Temp   = float(smash_conf["Modi"]["Box"]["Temperature"])
    dt     = float(smash_conf["General"]["Delta_Time"])
    total_time = float(smash_conf["General"]["End_Time"]) - t_start
    # Count elastic collisions
    coll_file = os.path.join(path,'collisions_binary.bin')
    res = count_elastic_scat(coll_file, t_start)
    N_events = res['ev_num']
    N_coll   = res['elast_coll']
    outfile.write('%i %i %.4f %.1f %.1f %i %i %.4f %.6f\n' %
           (N_coll, N_events, total_time, V, sigma, N_part, N_test, Temp, dt))
outfile.write('# %s\n' % res['smash_ver'])

outfile.close()

