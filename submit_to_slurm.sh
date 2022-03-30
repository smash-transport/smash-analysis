#!/bin/bash

if [[ "$#" -ne 5 ]] && [[ "$#" -ne 6 ]] && [[ "$#" -ne 7 ]]; then
    echo "Submit SLURM job for given target."
    echo
    echo "Usage: $0 SMASH_SRC_DIR EIGEN_SRC_DIR GSL_DIR PYTHIA_DIR TARGET [OUTPUT_DIR] [SAMPLED_LISTS]"
    exit 1
fi

path=$(readlink -f $0)
root=$(dirname $path)

target=$5
opt_argument=$6

if [[ "$#" -eq 5 ]]; then
    dt=`date +%Y_%m_%d_%H:%M:%S`
    output=$root/slurm_results_${dt}_${target}
elif [[ "$#" -eq 6 ]] && [[ "${opt_argument:0:4}" = "http" ]]; then
    dt=`date +%Y_%m_%d_%H:%M:%S`
    output=$root/slurm_results_${dt}_${target}
    sampled_lists=$6
else
    output=$6
    sampled_lists=$7
fi

mkdir -p $output
cd $output

if [ ${target} = "all" ]
then
  sbatch --job-name=smash-ang_dist $root/slurm-job.sh angular_distributions $1 $2 $3 $4 $root $output/angular_distributions
  sbatch --job-name=smash-xsec $root/slurm-job.sh cross_sections $1 $2 $3 $4 $root $output/cross_sections
  sbatch --job-name=smash-det_bal $root/slurm-job.sh detailed_balance $1 $2 $3 $4 $root $output/detailed_balance
  sbatch --job-name=smash-dil $root/slurm-job.sh dileptons $1 $2 $3 $4 $root $output/dileptons
  sbatch --job-name=smash-el_box $root/slurm-job.sh elastic_box $1 $2 $3 $4 $root $output/elastic_box
  sbatch --job-name=smash-e_scan $root/slurm-job.sh energy_scan $1 $2 $3 $4 $root $output/energy_scan
  sbatch --job-name=smash-FOPI $root/slurm-job.sh FOPI_pions $1 $2 $3 $4 $root $output/FOPI_pions
  sbatch --job-name=smash-dens $root/slurm-job.sh densities $1 $2 $3 $4 $root $output/densities
  sbatch --job-name=smash-aburner $root/slurm-job.sh afterburner $1 $2 $3 $4 $root $output/afterburner $sampled_lists
else
  sbatch --job-name=smash-$target $root/slurm-job.sh $target $1 $2 $3 $4 $root $output
fi
