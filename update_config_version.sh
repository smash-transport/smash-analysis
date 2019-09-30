#!/bin/sh
if [[ "$#" -ne 1 ]]; then
    echo "Update all config version numbers."
    echo
    echo "Usage: $0 SMASH_VERSION_NUMBER"
    exit 1
fi

smash_version=$1
find ./test/ -maxdepth 3 -name *.yaml | xargs sed -i '1 s/^.*$/Version: '$smash_version'/'

echo 'Updated all config version numbers successfully to' $smash_version
