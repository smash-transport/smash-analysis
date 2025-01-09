#!/bin/bash

if [[ "$#" -ne 2 ]] && [[ "$#" -ne 3 ]] && [[ "$#" -ne 4 ]]; then
    echo "Submit SLURM job for given target."
    echo
    echo "Usage: $0 CONTAINER TARGET [OUTPUT_DIR] [SAMPLED_LISTS]"
    exit 1
fi

path=$(readlink -f $0)
root=$(dirname $path)

target=$2
opt_argument=$3

if [[ "$#" -eq 2 ]]; then
    dt=`date +%Y_%m_%d_%H:%M:%S`
    output=$root/slurm_results_${dt}_${target}
elif [[ "$#" -eq 3 ]] && [[ "${opt_argument:0:4}" = "http" ]]; then
    dt=`date +%Y_%m_%d_%H:%M:%S`
    output=$root/slurm_results_${dt}_${target}
    sampled_lists=$3
else
    output=$3
    sampled_lists=$4
fi

mkdir -p $output
cd $output

if [ ${target} = "all" ]
then
  sbatch --job-name=smash-ang_dist $root/slurm-job-container.sh angular_distributions $1  $root $output/angular_distributions
  sbatch --job-name=smash-xsec $root/slurm-job-container.sh cross_sections $1  $root $output/cross_sections
  sbatch --job-name=smash-det_bal $root/slurm-job-container.sh detailed_balance $1  $root $output/detailed_balance
  sbatch --job-name=smash-dil $root/slurm-job-container.sh dileptons $1  $root $output/dileptons
  sbatch --job-name=smash-el_box $root/slurm-job-container.sh elastic_box $1  $root $output/elastic_box
  sbatch --job-name=smash-e_scan $root/slurm-job-container.sh energy_scan $1  $root $output/energy_scan
  sbatch --job-name=smash-FOPI $root/slurm-job-container.sh FOPI_pions $1  $root $output/FOPI_pions
  sbatch --job-name=smash-dens $root/slurm-job-container.sh densities $1  $root $output/densities
  sbatch --job-name=smash-aburner $root/slurm-job-container.sh afterburner $1  $root $output/afterburner $sampled_lists
else
  sbatch --job-name=smash-$target $root/slurm-job-container.sh $target $1 $root $output
fi
