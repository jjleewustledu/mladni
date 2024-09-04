#!/bin/bash

# Use absolute paths!
inFiles="/scratch/jjlee/Singularity/ADNI/NMF_FDG/warped_dementia_amyneg/nifti_files.csv"
outputDir="/scratch/jjlee/Singularity/ADNI/NMF_FDG/warped_dementia_amyneg"
scriptDir="/home/aris_data/ADNI_FDG/VolBin"
ranks="2 4 8" ###"2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40"

echo ${inFiles}
echo ${outputDir}

# submit master script

# number of components to evaluate
for b in $ranks
do
	mkdir -p ${outputDir}/NumBases${b}
	${scriptDir}/create_submit_extractBasesMT.sh ${inFiles} ${b} ${outputDir} ${scriptDir}> ${outputDir}/NumBases${b}/submit_extractBasesMT.sh
	cd ${outputDir}/NumBases${b}
	echo $(pwd)
	sbatch ./submit_extractBasesMT.sh
done	
