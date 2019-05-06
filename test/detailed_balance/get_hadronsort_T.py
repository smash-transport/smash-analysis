import sys
import numpy as np
import smash_basic_scripts as sb

arg = sys.argv
file_to_analyze = arg[1]
output_file = arg[2]
pdg_string = arg[3]
tstart = float(arg[4])

pdg_list = map(int, pdg_string.split(','))
N_pdgs = len(pdg_list)

with sb.BinaryReader(file_to_analyze) as reader:
    smash_version = reader.smash_version
    format_version = reader.format_version

    # Set histgram binnings
    Ebins = np.linspace(0.0, 5., 101)

    # Prepare empty lists of histograms
    Ehist = []
    for i in xrange(N_pdgs):
        h, bins = np.histogram([], bins = Ebins)
        Ehist.append(h)

    # Read file and fill histograms
    for block in reader:
        if (block['type'] == 'p'):
            time = sb.get_block_time(block)
            if (time > tstart):
                E   = block['part']['p'][:,0]
                pdg = block['part']['pdgid']
                for i in xrange(N_pdgs):
                    Ehist_ev, bins =  np.histogram(E, bins = Ebins,
                                        weights = (pdg == pdg_list[i]).astype(int))
                    Ehist[i] += Ehist_ev
        if (block['type'] == 'f'):  # end of event
            event_num = block['nevent'] + 1
            if (event_num == 750): break
            # print "finished analyzing event ", event_num
        if (block['type'] == 'i'):
            print "Error: interaction block in particles file."
            sys.exit(1)

# Get bin centers instead of edges
Ebins = 0.5*(Ebins[:-1] + Ebins[1:])

# Write output
with open(output_file,'w') as f:
    f.write("# Number of events\n")
    f.write("%d\n" % event_num)
    f.write("# Start counting at time [fm/c]\n")
    f.write("%.3f\n" % tstart)
    f.write("# Pdg list\n")
    f.write("%s\n" % pdg_string)

with open(output_file, 'a') as f:
    for i in xrange(N_pdgs):
        f.write("# Pdg %d: (E N(E))\n" % pdg_list[i])
        np.savetxt(f, Ebins, newline = ' ', fmt = '%.3f')
        f.write("\n")
        np.savetxt(f, Ehist[i], newline = ' ', fmt = '%.3f')
        f.write("\n")

