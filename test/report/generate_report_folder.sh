#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Expected 3 arguments, received $#."
    echo
    echo "Usage: bash $0 SMASH_VERSION PDF_CONFIG_DIR REPORT_DIR"
    exit 1
fi

smash_version=$1
collection_dir=$2       # contains all available PDFs and configs of the analysis
                        # suite run.
report_dir=$3           # report directory for html page


mkdir ${report_dir}
echo "Found results folder: ${collection_dir}"
echo "Will generate report folder: ${report_dir}"

#store location of generate_report_folder.sh file
current_dir=`pwd`

cd ${collection_dir}
for d in */ ; do
  test_name=`basename "$d"`
  cd ${current_dir}
  echo "Generating html file for: ${test_name}"
  python2 ./generate_html.py ${collection_dir}/$d --version $smash_version -o ${report_dir}/${test_name}
  cd ${collection_dir}
done

cd ${current_dir}
echo "Generating html file for: frontpage"
python2 generate_html.py ${collection_dir} --version $smash_version -o ${report_dir} --frontpage
