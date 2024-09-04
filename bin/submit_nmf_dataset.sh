#!/bin/bash

set -e
WORK=/scratch/jjlee/Singularity/ADNI/NMF_FDG
usage="$(basename "$0") [-h] [-d nmf_dataset]

where:  -h shows this help
        -d <nmf_dataset> specifies a subfolder of $WORK for NMF inference"

outputDir=$WORK
while getopts ':hd:' option; do
    case "$option" in
	h) echo "$usage"
	   exit
	   ;;
	d) outputDir=$WORK/$OPTARG
	   ;;
	:) printf "missing argument of -%s\n" "$OPTARG" >&2
	   echo "$usage" >&2
	   exit 1
	   ;;
       \?) printf "unknown option: -%s\n" "$OPTARG" >&2
	   echo "$usage" >&2
	   exit 1
	   ;;
    esac
done
shift $((OPTIND - 1))

if [ ! -d "$WORK" ]; then
    echo "ERROR:  $WORK not found" >&2
    exit 1    
fi

if [ ! -d "$outputDir" ]; then
    mkdir $outputDir
fi
echo "outputDir is $outputDir"

inFiles="$outputDir/nifti_files.csv"
if [ ! -f "$inFiles" ]; then
    echo "ERROR:  $inFiles not found" >&2
    exit 1
fi

#=================================================================================================

scriptDir="/home/aris_data/ADNI_FDG/VolBin"
ranks="2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40" # number of components to evaluate

# submit master script
for b in $ranks
do
	mkdir -p ${outputDir}/NumBases${b}
	${scriptDir}/create_submit_extractBasesMT.sh ${inFiles} ${b} ${outputDir} ${scriptDir} > ${outputDir}/NumBases${b}/submit_extractBasesMT.sh
	cd ${outputDir}/NumBases${b}
	echo "submitting slurm job from $(pwd)"
	sbatch ./submit_extractBasesMT.sh
done	
