#!/bin/bash
#SBATCH --mem-per-cpu=2000M
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --partition=long
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=<goetz@itp.uni-frankfurt.de>
#SBATCH --mail-type=ALL


if [[ "$#" -ne 4 ]] && [[ "$#" -ne 5 ]]; then
    echo "Expected 4 to 5 arguments, got $#."
    echo
    echo "Usage: sbatch $0 CMAKE_TARGET CONTAINER ANALYSIS_SRC_DIR OUTPUT_DIR [SAMPLED LISTS]"
    exit 1
fi

target=$1
container=$2
analysis_dir=$3
output_dir=$4
sampled_lists=$5


build_dir=$output_dir/build

echo "Running on ${SLURM_NNODES} nodes."
echo "Number of tasks: ${SLURM_NTASKS}"
echo "Number of CPUs per node: ${SLURM_CPUS_ON_NODE}"
echo "Executed from: ${SLURM_SUBMIT_DIR}"
echo "List of nodes: ${SLURM_JOB_NODELIST}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Analysis target: ${target}"
echo "Analysis directory: ${analysis_dir}"
echo "Output directory: ${output_dir}"
echo "Build directory: ${build_dir}"

date

singularity exec $container bash -c "mkdir -p $output_dir \
&& export PYTHONPATH=\"/lustre/hyihp/ngoetz/smash-analysis-extras/python_scripts/\"\
&& cd $output_dir \
&& cmake -DSMASH_PATH=\"/SMASH/smash_bin\" -B$output_dir -H$analysis_dir -DSAMPLED_LISTS=$sampled_lists -DEXP_DATA=\"/lustre/hyihp/ngoetz/smash-analysis-extras/experimental_data/\"\
&& \
if [ \"${target}\" = \"spectra\" ] || [ \"${target}\" = \"pp_collisions\" ] || [ \"${target}\" = \"FOPI_pions\" ]; then
  make ${target}_sims -j$SLURM_CPUS_ON_NODE \
  && make ${target}_plots -j$SLURM_CPUS_ON_NODE
elif [ \"${target}\" = \"elastic_box\" ] || [ \"${target}\" = \"detailed_balance\" ] || [ \"${target}\" = \"angular_distributions\" ] || [ \"${target}\" = \"cross_sections\" ] || [ \"${target}\" = \"energy_scan\" ] || [ \"${target}\" = \"afterburner\" ] || [ \"${target}\" = \"densities\" ]; then
  make ${target}_sims -j$SLURM_CPUS_ON_NODE \
  && make ${target}_analysis -j$SLURM_CPUS_ON_NODE \
  && make ${target}_plots
else
  make $target -j$SLURM_CPUS_ON_NODE
fi"


date
