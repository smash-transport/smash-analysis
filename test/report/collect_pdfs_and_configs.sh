#!/bin/bash

if [[ "$#" -ne 2 ]] && [[ "$#" -ne 3 ]]; then
    echo "Store .txt and .dat files for future comparison and PDFs for wiki page"
    echo
    echo "Usage: $0 SMASH_VERSION RESULTS_DIRECTORY [PDF_CONFIG_DIRECTORY]"
    exit 1
fi

smash_version=$1
results=$2     # e.g. /lustre/.../smash-analysis-devel

# directory where all PDFs and corresponding configs are stored
if [ "$#" -eq 3 ]; then
  collection_dir=$3
else
  collection_dir="${results}/PDFs_configs_${smash_version}"
fi

# prepare collection directory
if [ ! -d "$collection_dir" ]; then
    mkdir -p "$collection_dir"
else
    rm -rf -- "$collection_dir"
    mkdir -p "$collection_dir"
    echo "WARNING: Directory to collect PDFs and configs already existed, it is being overwritten."
fi

# create top-level category dirs (this is what you want to see in ls)
mkdir -p "${collection_dir}/FOPI_pions"
mkdir -p "${collection_dir}/afterburner"
mkdir -p "${collection_dir}/angular_distributions"
mkdir -p "${collection_dir}/cross_sections"
mkdir -p "${collection_dir}/densities"
mkdir -p "${collection_dir}/detailed_balance"
mkdir -p "${collection_dir}/dileptons"
mkdir -p "${collection_dir}/elastic_box"
mkdir -p "${collection_dir}/energy_scan"

echo ""
echo "Copying pdf files now ..."

########################
#        PDFs          #
########################

# angular_distributions
cd "${results}/angular_distributions/test/angular_distributions"
find . -maxdepth 3 -name "*.pdf" \
  | xargs -r cp --parents -t "${collection_dir}/angular_distributions"

# cross_sections
cd "${results}/cross_sections/test/cross_sections"
find . -maxdepth 3 -name "*.pdf" \
  | xargs -r cp --parents -t "${collection_dir}/cross_sections"

# detailed_balance
cd "${results}/detailed_balance/test/detailed_balance"
find . -maxdepth 3 -name "*.pdf" \
  | xargs -r cp --parents -t "${collection_dir}/detailed_balance"

# dileptons
cd "${results}/dileptons/test/dileptons"
find . -maxdepth 3 -name "plot_mass.pdf" \
  | xargs -r cp --parents -t "${collection_dir}/dileptons"

# elastic_box
cd "${results}/elastic_box/test/elastic_box"
find . -maxdepth 2 -name "*.pdf" \
  | xargs -r cp --parents -t "${collection_dir}/elastic_box"

# FOPI_pions
cd "${results}/FOPI_pions/test/FOPI_pions"
find . -maxdepth 2 -name "*.pdf" \
  | xargs -r cp --parents -t "${collection_dir}/FOPI_pions"

# ENERGY SCAN PDFs – reordered into subfolders
mkdir -p "${collection_dir}/energy_scan/mtspectra"
mkdir -p "${collection_dir}/energy_scan/yspectra"
mkdir -p "${collection_dir}/energy_scan/meanmt"
mkdir -p "${collection_dir}/energy_scan/meanpt"
mkdir -p "${collection_dir}/energy_scan/total_multiplicity"
mkdir -p "${collection_dir}/energy_scan/midrapidity_yield"

cd "${results}/energy_scan/test/energy_scan"
find . -maxdepth 2 -name "yspectra*.pdf" \
  | xargs -r cp -t "${collection_dir}/energy_scan/yspectra/"
find . -maxdepth 2 -name "mtspectra*.pdf" \
  | xargs -r cp -t "${collection_dir}/energy_scan/mtspectra/"
find . -maxdepth 2 -name "meanmt*.pdf" \
  | xargs -r cp -t "${collection_dir}/energy_scan/meanmt/"
find . -maxdepth 2 -name "meanpt*.pdf" \
  | xargs -r cp -t "${collection_dir}/energy_scan/meanpt/"
find . -maxdepth 2 -name "midrapidity_yield*.pdf" \
  | xargs -r cp -t "${collection_dir}/energy_scan/midrapidity_yield/"
find . -maxdepth 2 -name "total_multiplicity*.pdf" \
  | xargs -r cp -t "${collection_dir}/energy_scan/total_multiplicity/"

# AFTERBURNER PDFs – reordered into subfolders
mkdir -p "${collection_dir}/afterburner/mtspectra"
mkdir -p "${collection_dir}/afterburner/yspectra"
mkdir -p "${collection_dir}/afterburner/meanmt"
mkdir -p "${collection_dir}/afterburner/meanpt"
mkdir -p "${collection_dir}/afterburner/total_multiplicity"
mkdir -p "${collection_dir}/afterburner/midrapidity_yield"
mkdir -p "${collection_dir}/afterburner/LHC"
mkdir -p "${collection_dir}/afterburner/RHIC"

cd "${results}/afterburner/test/afterburner"
find . -maxdepth 2 -name "yspectra*.pdf" \
  | xargs -r cp -t "${collection_dir}/afterburner/yspectra/"
find . -maxdepth 2 -name "mtspectra*.pdf" \
  | xargs -r cp -t "${collection_dir}/afterburner/mtspectra/"
find . -maxdepth 2 -name "meanmt*.pdf" \
  | xargs -r cp -t "${collection_dir}/afterburner/meanmt/"
find . -maxdepth 2 -name "meanpt*.pdf" \
  | xargs -r cp -t "${collection_dir}/afterburner/meanpt/"
find . -maxdepth 2 -name "midrapidity_yield*.pdf" \
  | xargs -r cp -t "${collection_dir}/afterburner/midrapidity_yield/"
find . -maxdepth 2 -name "total_multiplicity*.pdf" \
  | xargs -r cp -t "${collection_dir}/afterburner/total_multiplicity/"

# densities
cd "${results}/densities/test/densities"
find . -maxdepth 2 -name "*.pdf" \
  | xargs -r cp --parents -t "${collection_dir}/densities"

echo "Done copying pdf files."
echo ""
echo "Copying config files now ..."

########################
#       CONFIGS        #
########################

# angular_distributions
cd "${results}/angular_distributions/test/angular_distributions"
find . -maxdepth 4 -name "config.yaml" \
  | xargs -r cp --parents -t "${collection_dir}/angular_distributions"

# cross_sections
cd "${results}/cross_sections/test/cross_sections"
find . -maxdepth 4 -name "config.yaml" \
  | xargs -r cp --parents -t "${collection_dir}/cross_sections"

# detailed_balance
cd "${results}/detailed_balance/test/detailed_balance"
find . -maxdepth 4 -name "config.yaml" \
  | xargs -r cp --parents -t "${collection_dir}/detailed_balance"

# dileptons
cd "${results}/dileptons/test/dileptons"
find . -maxdepth 4 -name "config.yaml" \
  | xargs -r cp --parents -t "${collection_dir}/dileptons"

# elastic_box
cd "${results}/elastic_box/test/elastic_box"
find . -maxdepth 4 -name "config.yaml" \
  | xargs -r cp --parents -t "${collection_dir}/elastic_box"

# FOPI_pions
cd "${results}/FOPI_pions/test/FOPI_pions"
find . -maxdepth 4 -name "config.yaml" \
  | xargs -r cp --parents -t "${collection_dir}/FOPI_pions"

# ENERGY SCAN configs
cd "${results}/energy_scan/test/energy_scan"
find . -maxdepth 5 -name "config.yaml" \
  | xargs -r cp --parents -t "${collection_dir}/energy_scan"

# AFTERBURNER configs (special LHC/RHIC handling)
cd "${results}/afterburner/test/afterburner"
find . -maxdepth 5 -path "*/LHC/1/data/config.yaml" \
  | xargs -r cp -t "${collection_dir}/afterburner/LHC"
find . -maxdepth 5 -path "*/RHIC/1/data/config.yaml" \
  | xargs -r cp -t "${collection_dir}/afterburner/RHIC"

# densities
cd "${results}/densities/test/densities"
find . -maxdepth 4 -name "config.yaml" \
  | xargs -r cp --parents -t "${collection_dir}/densities"

echo "Done copying config files."
echo ""
echo "Finished copying PDFs and configs."
echo "Collected into: ${collection_dir}"

