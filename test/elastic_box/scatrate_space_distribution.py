import os
import sys
import yaml
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../../python_scripts')
import smash_basic_scripts as sb
import numpy as np

smash_input = sys.argv[1]
yaml_input  = sys.argv[2]
output_file = sys.argv[3]

bfile = open(smash_input, 'rb')
smash_version, format_version = sb.read_binary_header(bfile)
with open(yaml_input, 'r') as cf:
    smash_conf = yaml.load(cf, Loader=yaml.loader.BaseLoader)
parameters = dict()
L = float(smash_conf["Modi"]["Box"]["Length"])
parameters['$V$ [fm$^3$]'] = str(L*L*L)
parameters['N'] = smash_conf["Modi"]["Box"]["Init_Multiplicities"]["111"]
parameters['$\sigma$ [mb]'] = smash_conf["Collision_Term"]["Sigma"]
parameters['$N_{test}$'] = smash_conf["General"]["Testparticles"]
parameters['T [GeV]'] = smash_conf["Modi"]["Box"]["Temperature"]
parameters['dt [fm/c]'] = smash_conf["General"]["Delta_Time"]
parameters['time [fm/c]'] = smash_conf["General"]["End_Time"]
total_time = float(parameters['time [fm/c]'])

# Lists of collision coordinates and times
t = []
x = []
y = []
z = []
while True:
    block = sb.read_binary_block(bfile)
    if (block == None):  # end of file
        break
    if (block['type'] == 'f'):  # end of event
        event_num = block['nevent'] + 1
        if (event_num % 20 == 0):
            print "finished analyzing event ", event_num
    if (block['type'] == 'i' and sb.is_elastic22(block)): # 2->2 elastic
            coll_time = sb.get_block_time(block)
            #if (block['incoming']['r'][:,1].max() > L):
            #    print "x out-of-box: ", block['incoming']['r'][:,1].max()
            #if (block['incoming']['r'][:,1].min() < 0.0):
            #    print "x out-of-box: ", block['incoming']['r'][:,1].min()
            coll_x = block['incoming']['r'][:,1].mean()
            coll_y = block['incoming']['r'][:,2].mean()
            coll_z = block['incoming']['r'][:,3].mean()
            if (coll_x < 0.0): coll_x +=L
            if (coll_x > L): coll_x -=L
            if (coll_y < 0.0): coll_y +=L
            if (coll_y > L): coll_y -=L
            if (coll_z < 0.0): coll_z +=L
            if (coll_z > L): coll_z -=L
            t.append(coll_time)
            x.append(coll_x)
            y.append(coll_y)
            z.append(coll_z)
bfile.close()
parameters['$N_{ev}$'] = str(event_num)
parameters['version'] = smash_version[0]

# Create histograms: number of collisions binned by time and x
t_hist, t_bin_edges = np.histogram(t, bins = 100, range = (0, total_time))
x_hist, x_bin_edges = np.histogram(x, bins = 20, range = (0, L))
y_hist, y_bin_edges = np.histogram(y, bins = 20, range = (0, L))
z_hist, z_bin_edges = np.histogram(z, bins = 20, range = (0, L))
t_bins = 0.5 * (t_bin_edges[:-1] + t_bin_edges[1:])
x_bins = 0.5 * (x_bin_edges[:-1] + x_bin_edges[1:])
y_bins = 0.5 * (y_bin_edges[:-1] + y_bin_edges[1:])
z_bins = 0.5 * (z_bin_edges[:-1] + z_bin_edges[1:])

# write results
with open(output_file, 'w') as f:
    f.write('#')
    yaml.safe_dump(parameters, f, width = 1000)
    np.savetxt(f, (t_bins, t_hist), fmt = '%.3f')
    np.savetxt(f, (x_bins, x_hist), fmt = '%.3f')
    np.savetxt(f, (y_bins, y_hist), fmt = '%.3f')
    np.savetxt(f, (z_bins, z_hist), fmt = '%.3f')
