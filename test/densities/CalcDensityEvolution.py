#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
sys.path.append(
    os.path.dirname(
        os.path.abspath(__file__)) +
    '/../../python_scripts')
from read_binary_output import *
import smash_basic_scripts as sb
import yaml

class DensityEvolution:
    def __init__(self, Path, NDim, NFolder, BoxLength, TimeStep):
        self.Path = Path
        self.NDim = NDim
        self.NFolder = NFolder
        self.BoxLength = BoxLength
        self.TimeStep = TimeStep

    def initialize_grid(self):
        grid = np.zeros((self.NDim, self.NDim), dtype=float)
        return grid

    def find_grid_index(self, x, y):
        indices = []
        for i in range(len(x)):
            x_index = int(x[i] * self.NDim / self.BoxLength)
            y_index = int(y[i] * self.NDim / self.BoxLength)
            indices.append((x_index, y_index))
        return indices

    def update_grid(self):
        Grids = []
        for evaluation_time in self.TimeStep:
            density_grid = self.initialize_grid()
            event_counter = 0
            for i in range(self.NFolder + 1):
                PFile = self.Path + str(i) + '/particles_binary.bin'
                with BinaryReader(PFile) as reader:
                    smash_version = reader.smash_version
                    block = reader.read_block()
                    while block is not None:
                        t = block['type']
                        assert(t == 'p' or t == 'f')
                        if (t == 'p'):
                            timestep = block['part']['r'][:,0][0]
                            if timestep == evaluation_time:
                                event_counter += 1
                                x = block['part']['r'][:,1]
                                y = block['part']['r'][:,2]
                                z = block['part']['r'][:,3]
                                particle_indices = self.find_grid_index(x, y)
                                for tpl in particle_indices:
                                    density_grid[tpl] += 1
                        block = reader.read_block()
            density_grid /= (event_counter)
            density_grid = density_grid.reshape(self.NDim**2)
            Grids.append(density_grid)
        return Grids, smash_version
 ##############################################################################
# with open(args.config_file, 'r') as f: d = yaml.load(f)


Loc = sys.argv[1]
config = sys.argv[2]
OutFile = sys.argv[3]
# DEFAULT Values
dim_grid = 40
NFolders = 10
with open(config, 'r') as f: config_file = yaml.safe_load(f)
TimeSteps = config_file['Output']['Output_Interval']
FinalTime = config_file['General']['End_Time']
Times = np.arange(0, FinalTime + TimeSteps, TimeSteps)
BoxLength = config_file['Modi']['Box']['Length']
RunClass = DensityEvolution(Loc, dim_grid, NFolders, BoxLength, Times)
density_grid, smash_version = RunClass.update_grid()
smash_analysis_version = sb.analysis_version_string()

write_out = open(OutFile, 'w')
write_out.write('# smash and smash analysis version\n')
write_out.write('{} {}\n'.format(smash_version, smash_analysis_version))
write_out.write('# Box length\n')
write_out.write('{}\n'.format(BoxLength))
write_out.write('# Number of grid cells along one direction\n')
write_out.write('{}\n'.format(dim_grid))
write_out.write('# Time steps\n')
for t in Times:
    write_out.write('{}'.format(t) + '\t')
write_out.write('\n')
write_out.write('# Grids\n')
for n in range(dim_grid**2):
    for i in range(len(Times)):
        write_out.write(str(density_grid[i][n]) + '\t')
    write_out.write('\n')
write_out.close()
