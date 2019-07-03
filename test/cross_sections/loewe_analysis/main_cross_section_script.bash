#!/bin/bash

# folder containing smash_basic_scripts.py
smashpython="../../python_scripts"
# folder containing SMASH executable
smashbuild="/path/to/SMASH/build"
# folder for running simulations and storing data
workpath="/path/to/data"
# path to cross section script
xsectionscript="../geometrical_cross_section.py"

nevents=50000

# get absolute paths
smashpython=$(readlink -f ${smashpython})
xsectionscript=$(readlink -f ${xsectionscript})
# add path to python scripts to includes known to python
export PYTHONPATH=$PYTHONPATH:$smashpython

# find the path to python for this system
pythonpath=$(which python2)
# update the python path in scripts
# replace only the first occurrence of shebang
# since we're dealing with paths with '/'s, let's use @ as delimiter for sed
sed -i "0,/#!/{s@#!.*@#!${pythonpath}@g}" ${xsectionscript}

# function for sending analysis jobs to SLURM
sbatchjob () {
    # Input parameters:
    # 1: Initial state identifier string
    # 2: PDG code of the first initial particle
    # 3: PDG code of the second initial particle
    # 4: Minimum collision energy sqrt(s) (in GeV)
    # 5: Maximum collision energy sqrt(s) (in GeV)
    # 6: Width of the energy bins (in GeV)
    local datafolder="${workpath}/${1}"
    # Overwrite existing data? 0 = False, 1 = True
    overwrite=1
    sbatch --job-name=xs_${1} --output=${1}_out.log --error=${1}_err.log simulate_collider.bash ${smashbuild} ${datafolder} ${xsectionscript} $2 $3 $4 $5 $6 $nevents $overwrite
    return 0
}

# Calculate sigma(sqrts)
# pi+ pi-
sbatchjob "pipi" 211 -211 0.3 1.2 0.05
# pp
sbatchjob "pp" 2212 2212 1.9 2.8 0.05
# np
sbatchjob "np" 2112 2212 1.9 2.8 0.05
# ppi+
sbatchjob "ppi+" 2212 211 1.1 1.6 0.02
# ppi-
sbatchjob "ppi-" 2212 -211 1.1 1.6 0.02
