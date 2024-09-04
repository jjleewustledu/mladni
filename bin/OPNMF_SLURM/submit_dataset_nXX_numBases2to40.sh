#!/bin/bash

inFiles="path to .csv file"
outputDir="OutputDirectory/"

echo ${inFiles}
echo ${outputDir}

scriptDir="path to VolBin"

# submit master script

# number of components to evaluate
for b in 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40
#for b in {1..20}
do
	mkdir -p ${outputDir}/NumBases${b}
	${scriptDir}/create_submit_extractBasesMT.sh ${inFiles} ${b} ${outputDir} > ${outputDir}/NumBases${b}/submit_extractBasesMT.sh
	cd ${outputDir}/NumBases${b}
	echo $(pwd)
	sbatch ./submit_extractBasesMT.sh
done	
