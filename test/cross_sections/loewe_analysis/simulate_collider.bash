#!/bin/bash

smashscript="smash_basic_scripts.py"
homepath="${PWD}"

# check the input arguments
if [ "$#" -lt 9 ]; then
    echo "ERROR: Expected 9 arguments, got $#"
    echo "Expected arguments:"
    echo "1. path to folder containing executable"
    echo "2. path where the simulations should be run"
    echo "3. name of analysis script"
    echo "4. PDG code of first particle"
    echo "5. PDG code of second particle"
    echo "6. lower limit of energy range"
    echo "7. upper limit of energy range"
    echo "8. energy step size"
    echo "9. number of events"
    exit 3
fi

execpath=$1
execpath=$(readlink -f $execpath)
if [ ! -f "${execpath}/smash" ]; then
    echo "ERROR: Did not find executable ${execpath}/smash"
    exit 3
fi

runpath=$2
runpath=$(readlink -m $runpath)
if mkdir -p ${runpath}; then
    echo "INFO: Running simulations in ${runpath}"
else
    echo "ERROR: Could not create ${runpath}"
    exit 3
fi

xsectionscript=$3
if [ ! -f "${xsectionscript}" ]; then
    echo "ERROR: Did not find required file ${xsectionscript}"
    exit 3
fi
pdg1=$4
pdg2=$5
emin=$6
emax=$7
ebin=$8
nevents=$9

overwrite=0
if [ "$#" -eq 10 ]; then
    overwrite=${10}
fi

# inform user if existing data is overwritten or not
if [ "$overwrite" -eq 0 ]; then
    echo "INFO: Existing data in ${runpath} will not be overwritten."
    echo "Simulations will be skipped in such cases."
else
    echo "INFO: Existing data in ${runpath} is overwritten."
fi

# prepare run directory for the simulations
cp ${execpath}/smash ${runpath}
cp ${execpath}/particles.txt ${runpath}
cp ${execpath}/decaymodes.txt ${runpath}
cp cross_section_config.yaml ${runpath}/config.yaml
cp ${xsectionscript} ${runpath}

energyvalues=()
nbins=$(echo "(${emax}-${emin})/${ebin}" | bc)
nevalues=0
while [ "${nevalues}" -le "$nbins" ]
do
    energyvalues+=($(echo "${nevalues}*${ebin}+${emin}" | bc -l))
    ((nevalues++))
done

cd ${runpath}
sed -i -e "s/EVENTVALUE/${nevents}/g" config.yaml
sed -i -e "s/PROJECTILEPDG/${pdg1}/g" config.yaml
sed -i -e "s/TARGETPDG/${pdg2}/g" config.yaml

#SBATCH --partition=parallel
#SBATCH --constraint=intel20
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=00:55:00
#SBATCH --array=0-19:20

for I in {0..19}
do
    if [ "$I" -lt "${nevalues}" ]; then
	E=${energyvalues[I]}
	rseed=$I
	if [ ! -f "./sqrts_${E}/collisions_binary.bin" ] || [ "$overwrite" -ne 0 ]; then
	    time -p ./smash -o "sqrts_${E}" -c "General: {Randomseed: ${rseed}}" -c "Modi: {Nucleus: {Sqrtsnn: ${E}}}" > smash_${E}.out &
	fi
    fi
done

wait

if [ "${nevalues}" -gt 20 ]; then
    for I in {20..39}
    do
	if [ "$I" -lt "${nevalues}" ]; then
	    E=${energyvalues[I]}
	    rseed=$I
	    if [ ! -f "./sqrts_${E}/collisions_binary.bin" ] || [ "$overwrite" -ne 0 ]; then
		time -p ./smash -o "sqrts_${E}" -c "General: {Randomseed: ${rseed}}" -c "Modi: {Nucleus: {Sqrtsnn: ${E}}}" > smash_${E}.out &
	    fi
	fi
    done
fi
    
wait

scriptname=$(basename $xsectionscript)
# calculate the cross sections and write the results to file
resultsfile="xsections_${pdg1}_${pdg2}_ss_${emin}_${emax}.data"
echo "INFO: Simulations done, running analysis."
echo "# Analysis done with ${scriptname}" > $resultsfile
time -p ./${scriptname} -f -b $ebin $pdg1 $pdg2 ${runpath}/sqrts_*/*.bin >> $resultsfile
cp ${resultsfile} ${homepath}
cd ${homepath}
