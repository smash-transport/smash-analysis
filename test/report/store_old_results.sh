#!/bin/bash

if [[ "$#" -ne 2 ]]; then
    echo "Store txt and dat files for future comparison and PDFs for wiki page"
    echo
    echo "Usage: $0 SMASH_VERSION RESULTS_DIRECTORY"
    exit 1
fi

smash_version=$1
results=$2             # directory where results can be found

store_dir="${results}/old_results_${smash_version}"  # for future comparison

if [ ! -d $store_dir ]; then
    mkdir ${store_dir}
else
    echo "Directory for old results already exists, please remove or rename."
    exit
fi

echo "Copying old data files now ..."
# copy only .txt and .dat files that are relevant for future comprison
# folder structure is preserved
cd ${results}/angular_distributions/test
find angular_distributions -maxdepth 3 -name "*.dat" | xargs cp --parents -t ${store_dir}/
cd ${results}/cross_sections/test
find cross_sections -maxdepth 3 -name "sqrts_totalXSec.txt" | xargs cp --parents -t ${store_dir}/
cd ${results}/detailed_balance/test
find detailed_balance -maxdepth 3 -name "Nreact_by_Nisopingroup.txt" | xargs cp --parents -t ${store_dir}/
cd ${results}/dileptons/test
find dileptons -maxdepth 3 -name "dN_dm_tot.txt" | xargs cp --parents -t ${store_dir}/
mkdir ${store_dir}/elastic_box
cd ${results}/elastic_box/test
find elastic_box -maxdepth 2 -name "scatrate_vs_*.txt" | xargs cp -t ${store_dir}/elastic_box/
cd ${results}/FOPI_pions/test
find FOPI_pions -maxdepth 2 -name "*.txt" | xargs cp --parents -t ${store_dir}/
cd ${results}/energy_scan/test
find energy_scan -maxdepth 2 -name "*.txt" | xargs cp --parents -t ${store_dir}/
cd ${results}/densities/test
find densities -maxdepth 3 -name "*.dat" | xargs cp --parents -t ${store_dir}/


echo "Done copying old data files"
echo "Finished successfully."
