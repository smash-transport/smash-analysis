#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --mem=MaxMemPerNode
#SBATCH --partition=long

if [ "$#" -ne 6 ]; then
    echo "Expected 6 arguments, got $#."
    echo
    echo "Usage: sbatch $0 CMAKE_TARGET SMASH_SRC_DIR EIGEN_SRC_DIR PYTHIA_DIR ANALYSIS_SRC_DIR OUTPUT_DIR"
    exit 1
fi

target=$1
smash_dir=$2
eigen_dir=$3
pythia_dir=$4
analysis_dir=$5
output_dir=$6

#build_dir=$(mktemp -d /tmp/smash-build.XXXXXXXXXX)
build_dir=$output_dir/build

echo "Running on ${SLURM_NNODES} nodes."
echo "Number of tasks: ${SLURM_NTASKS}"
echo "Number of CPUs per node: ${SLURM_CPUS_ON_NODE}"
echo "Executed from: ${SLURM_SUBMIT_DIR}"
echo "List of nodes: ${SLURM_JOB_NODELIST}"
echo "Job id: ${SLURM_JOB_ID}"
echo "SMASH directory: ${smash_dir}"
echo "Analysis target: ${target}"
echo "Analysis directory: ${analysis_dir}"
echo "Output directory: ${output_dir}"
echo "Build directory: ${build_dir}"
echo

date

mkdir -p $output_dir \
&& cd $smash_dir \
&& cmake -DCMAKE_BUILD_TYPE=Release -DPythia_CONFIG_EXECUTABLE=$pythia_dir/bin/pythia8-config -DCMAKE_INSTALL_PREFIX=$eigen_dir -B$build_dir -H$smash_dir \
&& cd $build_dir \
&& make smash -j$SLURM_CPUS_ON_NODE \
&& cd $output_dir \
&& cmake -DSMASH_PATH=$build_dir -B$output_dir -H$analysis_dir \
&&
if [ ${target} = "spectra" ] || [ ${target} = "pp_collisions" ] || [ ${target} = "FOPI_pions" ]
then
  make ${target}_sims -j$SLURM_CPUS_ON_NODE \
  && make ${target}_plots -j$SLURM_CPUS_ON_NODE
elif [ ${target} = "elastic_box" ] || [ ${target} = "detailed_balance" ] || [ ${target} = "angular_distributions" ] || [ ${target} = "cross_sections" ] || [ ${target} = "energy_scan" ]
then
  make ${target}_sims -j$SLURM_CPUS_ON_NODE \
  && make ${target}_analysis -j$SLURM_CPUS_ON_NODE \
  && make ${target}_plots
else
  make $target -j$SLURM_CPUS_ON_NODE
fi

date
