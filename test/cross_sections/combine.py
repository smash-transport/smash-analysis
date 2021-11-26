"""
Combine analysis output files for different energies into one file.

Note that it is important the there is exactly one file per energy.
"""

import sys
import os
import argparse

sys.path.append(
    os.path.dirname(
        os.path.abspath(__file__)) +
    '/../../python_scripts')
#from ordered_default_dict import OrderedDefaultDict, OrderedSet
from ordered_default_dict import OrderedSet
from collections import defaultdict

#reload(sys)
#sys.setdefaultencoding('utf-8')

def parse_arguments():
    """Parse and return the command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("output_file", help="Output of combined data")
    parser.add_argument("--data", nargs="+", help="Input data to be combined")
    args = parser.parse_args()
    return args

def store_results(results_file, x, y, version):
    store_file = open(results_file, "w")
    store_file.write("sqrts" + "\t" + "x_sec" + "\n")
    for i in range(0,len(x)):
        store_file.write(str(x[i]) + "\t" + str(y[i]) + "\n")
    store_file.write("# " + str(version))

def flatten(iterable):
    return [i for subiter in iterable for i in subiter]

if __name__ == '__main__':
    args = parse_arguments()
    cross_sections = defaultdict(float)
    print("In combinbe.py - main:")
    print(str(cross_sections))
    energies = OrderedSet()
    colnames = OrderedSet()

    with open(args.output_file, 'w') as out:
        first = True
        for d in args.data:
            with open(d, 'r') as f:
                if first:
                    # Read and write the header (should be the same for all files)
                    out.write(f.readline())
                    out.write(f.readline())
                    # extract SMASH version to store for old_results
                    header_line3 = f.readline()
                    version = header_line3.split()[1]
                    # continue writing header
                    out.write(header_line3)
                    first = False
                else:
                    # Discard the header (should be the same for all files)
                    f.readline()
                    f.readline()
                    f.readline()

                # We have to parse the columns, because they will be different
                # for each energy
                split_line = f.readline().rstrip().split(' ')
                start = split_line[0]
                current_colnames = split_line[1:]
                # Add the energy at the very beginning
                colnames.add(current_colnames[0])

                while True:
                    line = f.readline().rstrip()
                    if line:
                        values = line.split(' ')
                        energy = values[0]
                        energies.add(energy)
                        for i in range(1, len(values)):
                            colname = current_colnames[i]
                            if (energy, colname) in cross_sections:
                                print("In cycle: "+energy+"  "+colname+"  "+str(cross_sections[(energy, colname)]))
                                print('WARN: duplicate (energy, colname) = ({}, {}) when combining data\nold: {}, new: {}'.format(energy, colname, cross_sections[(energy, colname)], values[i]))
                                value = max(cross_sections[(energy, colname)], values[i])
                            else:
                                value = values[i]
                            colnames.add(colname)

                            cross_sections[(energy, colname)] = value
                    else:
                        break

        # Write cross section to output file
        energies = sorted(energies)
        colnames = list(colnames)
        old_colnames = colnames[:]
        # Sort everything except for energy and total and elastic if present
        # Make sure intermediate-state data is plotted last
        energy = colnames[0]
        assert 'sqrt' in energy
        # Temporarily get rid of error labels for column reordering
        old_len = len(colnames)
        colnames = [c for c in colnames if not c.endswith('_err')]
        assert float(len(old_colnames) - 1) / 2 + 1 == float(len(colnames))
        colnames.remove('total')
        if 'elastic' in colnames:
            colnames.remove('elastic')
            elastic = ['elastic']
        else:
            elastic = []
        with_arrow = []
        without_arrow = []
        for name in colnames[1:]:
            if '→' in name:
                with_arrow.append(name)
            else:
                without_arrow.append(name)
        # Reorder column names
        colnames = ['total'] + elastic + \
                   sorted(without_arrow) + sorted(with_arrow)
        # Add error labels again
        colnames_with_errors = flatten([(i, i + '_err') for i in colnames])
        colnames = [energy] + colnames_with_errors
        assert sorted(colnames) == sorted(old_colnames)
        out.write(start)
        out.write(' ')
        out.write(' '.join(colnames))
        out.write('\n')

        for energy in energies:
            print(energy, end=' ', file=out)
            # Output the cross sections (and skip energy)
            for colname in colnames[1:]:
                xs = cross_sections[(energy, colname)]
                print(xs, end=' ', file=out)
            print(file=out)

        # generate output file for comparison to previous version
        # just necessary for 1 of the 4 runs, randomly pick xs_final

        filename = os.path.basename(args.output_file)
        total_xs = []
        for energy in energies:
            total_x_sec = total_xs.append(cross_sections[(energy, 'total')])

        if filename == "xs_final_grouped.dat":
            store_results("sqrts_totalXSec.txt", energies, total_xs, version)
