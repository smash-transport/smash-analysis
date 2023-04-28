#!/bin/bash

if [[ "$#" -ne 1 ]] && [[ "$#" -ne 2 ]]; then
    echo "Create html-report for SMASH-analysis suite results."
    echo
    echo "Usage: $0 ANALYSIS_RESULTS_DIR [REPORT_DIR]"
    exit 1
fi

# directory in where all files of the SMASH-analysis runs are located
results_dir=$1

# determine SMASH version from (arbitrary chosen) config from angular distributions run
smash_version=$(python3 ./python_scripts/version.py)

# directory where the report will be generated (html page)
if [ "$#" -eq 2 ]; then
    report_dir=$2
else
    report_dir="${results_dir}/report_${smash_version}"
fi

# directory where all PDFs and corresponding configs are stored and sorted
collection_dir="${results_dir}/PDFs_configs_${smash_version}"

cd test/report
./collect_pdfs_and_configs.sh $smash_version $results_dir $collection_dir
wait
echo "--------"
./generate_report_folder.sh $smash_version $collection_dir $report_dir
echo "Done generating report."
