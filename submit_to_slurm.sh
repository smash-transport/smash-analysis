#!/bin/bash

if [[ "$#" -ne 4 ]] && [[ "$#" -ne 5 ]]; then
    echo "Submit SLURM job for given target."
    echo
    echo "Usage: $0 SMASH_SRC_DIR EIGEN_SRC_DIR PYTHIA_DIR TARGET [OUTPUT_DIR]"
    exit 1
fi

path=$(readlink -f $0)
root=$(dirname $path)

target=$4

if [ "$#" -eq 4 ]; then
    dt=`date +%Y_%m_%d_%H:%M:%S`
    output=$root/slurm_results_${dt}_${target}
else
    output=$5
fi

mkdir -p $output
cd $output

if [ ${target} = "all" ]
then
  sbatch --job-name=smash-ang_dist $root/slurm-job.sh angular_distributions $1 $2 $3 $root $output/angular_distributions
  sbatch --job-name=smash-xsec $root/slurm-job.sh cross_sections $1 $2 $3 $root $output/cross_sections
  sbatch --job-name=smash-det_bal $root/slurm-job.sh detailed_balance $1 $2 $3 $root $output/detailed_balance
  sbatch --job-name=smash-dil $root/slurm-job.sh dileptons $1 $2 $3 $root $output/dileptons
  sbatch --job-name=smash-el_box $root/slurm-job.sh elastic_box $1 $2 $3 $root $output/elastic_box
  sbatch --job-name=smash-e_scan $root/slurm-job.sh energy_scan $1 $2 $3 $root $output/energy_scan
  sbatch --job-name=smash-FOPI $root/slurm-job.sh FOPI_pions $1 $2 $3 $root $output/FOPI_pions
else
  sbatch --job-name=smash-$target $root/slurm-job.sh $target $1 $2 $3 $root $output
fi
