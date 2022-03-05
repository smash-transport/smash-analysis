#!/bin/bash

if [[ "$#" -ne 2 ]] && [[ "$#" -ne 3 ]]; then
    echo "Store .txt and .dat files for future comparison and PDFs for wiki page"
    echo
    echo "Usage: $0 SMASH_VERSION RESULTS_DIRECTORY [PDF_CONFIG_DIRECTORY]"
    exit 1
fi

smash_version=$1

# directory in where all files of the SMASH-analysis runs are located
results=$2

# directory where all PDFs and corresponding configs are stored
if [ "$#" -eq 3 ]; then
  collection_dir=$3
else
  collection_dir="${results}/PDFs_configs_${smash_version}"
fi

if [ ! -d $collection_dir ]; then
    mkdir ${collection_dir}
else
    rm -r ${collection_dir}
    mkdir ${collection_dir}
    echo "WARNING: Directory to collect PDFs and configs already existed, it is being overwritten."
fi

echo ""
echo "Copying pdf files now ..."

# create new subdirectories for energy scan (files are being reordered):
mkdir ${collection_dir}/energy_scan
mkdir ${collection_dir}/energy_scan/mtspectra
mkdir ${collection_dir}/energy_scan/yspectra
mkdir ${collection_dir}/energy_scan/meanmt
mkdir ${collection_dir}/energy_scan/meanpt
mkdir ${collection_dir}/energy_scan/total_multiplicity
mkdir ${collection_dir}/energy_scan/midrapidity_yield

# create new subdirectories for afterburner run (files are being reordered):
mkdir ${collection_dir}/afterburner
mkdir ${collection_dir}/afterburner/mtspectra
mkdir ${collection_dir}/afterburner/yspectra
mkdir ${collection_dir}/afterburner/meanmt
mkdir ${collection_dir}/afterburner/meanpt
mkdir ${collection_dir}/afterburner/total_multiplicity
mkdir ${collection_dir}/afterburner/midrapidity_yield

# copy all PDF files necessary to generate the wiki page
# folder structure is preserved (except for energy scan, where we reorder)
cd ${results}/angular_distributions/test
find angular_distributions -maxdepth 3 -name "*.pdf" | xargs cp --parents -t ${collection_dir}/
cd ${results}/cross_sections/test
find cross_sections -maxdepth 3 -name "*.pdf" | xargs cp --parents -t ${collection_dir}/
cd ${results}/detailed_balance/test
find detailed_balance -maxdepth 3 -name "*.pdf" | xargs cp --parents -t ${collection_dir}/
cd ${results}/dileptons/test
find dileptons -maxdepth 3 -name "plot_mass.pdf" | xargs cp --parents -t ${collection_dir}/
cd ${results}/elastic_box/test
find elastic_box -maxdepth 2 -name "*.pdf" | xargs cp --parents -t ${collection_dir}/
cd ${results}/FOPI_pions/test
find FOPI_pions -maxdepth 2 -name "*.pdf" | xargs cp --parents -t ${collection_dir}/
cd ${results}/energy_scan/test
find energy_scan -maxdepth 2 -name "yspectra*.pdf" | xargs cp -t ${collection_dir}/energy_scan/yspectra/
find energy_scan -maxdepth 2 -name "mtspectra*.pdf" | xargs cp -t ${collection_dir}/energy_scan/mtspectra/
find energy_scan -maxdepth 2 -name "meanmt*.pdf" | xargs cp -t ${collection_dir}/energy_scan/meanmt/
find energy_scan -maxdepth 2 -name "meanpt*.pdf" | xargs cp -t ${collection_dir}/energy_scan/meanpt/
find energy_scan -maxdepth 2 -name "midrapidity_yield*.pdf" | xargs cp -t ${collection_dir}/energy_scan/midrapidity_yield/
find energy_scan -maxdepth 2 -name "total_multiplicity*.pdf" | xargs cp -t ${collection_dir}/energy_scan/total_multiplicity/
cd ${results}/afterburner/test
find afterburner -maxdepth 2 -name "yspectra*.pdf" | xargs cp -t ${collection_dir}/afterburner/yspectra/
find afterburner -maxdepth 2 -name "mtspectra*.pdf" | xargs cp -t ${collection_dir}/afterburner/mtspectra/
find afterburner -maxdepth 2 -name "meanmt*.pdf" | xargs cp -t ${collection_dir}/afterburner/meanmt/
find afterburner -maxdepth 2 -name "meanpt*.pdf" | xargs cp -t ${collection_dir}/afterburner/meanpt/
find afterburner -maxdepth 2 -name "midrapidity_yield*.pdf" | xargs cp -t ${collection_dir}/afterburner/midrapidity_yield/
find afterburner -maxdepth 2 -name "total_multiplicity*.pdf" | xargs cp -t ${collection_dir}/afterburner/total_multiplicity/
cd ${results}/densities/test
find densities -maxdepth 2 -name "*.pdf" | xargs cp --parents -t ${collection_dir}/


echo "Done copying pdf files."
echo ""
echo "Copying config files now ..."
# copy all config files necessary to generate the wiki page
# folder structure is preserved (except for energy scan, see above)
cd ${results}/angular_distributions/test
find angular_distributions -maxdepth 4 -name "config.yaml" | xargs cp --parents -t ${collection_dir}/
cd ${results}/cross_sections/test
find cross_sections -maxdepth 4 -name "config.yaml" | xargs cp --parents -t ${collection_dir}/
cd ${results}/detailed_balance/test
find detailed_balance -maxdepth 4 -name "config.yaml" | xargs cp --parents -t ${collection_dir}/
cd ${results}/dileptons/test
find dileptons -maxdepth 4 -name "config.yaml" | xargs cp --parents -t ${collection_dir}/
cd ${results}/elastic_box/test
find elastic_box -maxdepth 4 -name "config.yaml" | xargs cp --parents -t ${collection_dir}/
cd ${results}/FOPI_pions/test
find FOPI_pions -maxdepth 4 -name "config.yaml" | xargs cp --parents -t ${collection_dir}/
cd ${results}/energy_scan/test
find energy_scan -maxdepth 5 -name "config.yaml" | xargs cp --parents -t ${collection_dir}/
cd ${results}/afterburner/test
find afterburner -maxdepth 5 -name "config.yaml" | xargs cp --parents -t ${collection_dir}/
cd ${results}/densities/test
find densities -maxdepth 4 -name "config.yaml" | xargs cp --parents -t ${collection_dir}/
echo "Done copying config files."
echo ""
echo "Finished copying PDFs and configs."
