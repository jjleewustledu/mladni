#!/bin/bash  
### Job Name  
#SBATCH -J extractBasesMT  
### Save output and error files in the following folder  
#SBATCH -o /scratch/jjlee/Singularity/ADNI/NMF_FDG/NumBases2/slurm_job_output/output.log  
#SBATCH -e /scratch/jjlee/Singularity/ADNI/NMF_FDG/NumBases2/slurm_job_output/error.log  
### Start job with priority - default is zero  
#SBATCH --priority=0  
### Select 1 nodes with 4 CPU  
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=4  
### Select 24Gb RAM  
#SBATCH --mem=128G  
### Select wall time  
#SBATCH -t 24:00:00  
  
# set up useful variables  
  
MCRROOT=/export/matlab/MCR/R2018b/v95  
BinDir=/home/aris_data/ADNI_FDG/VolBin 
extractBases=/home/aris_data/ADNI_FDG/VolBin/run_extractBasesMT.sh
  
# parse input variables  
inFiles=/scratch/jjlee/Singularity/ADNI/bids/derivatives/globbed_dementia_amyneg.csv	  
numBases=2  
outputDir=/scratch/jjlee/Singularity/ADNI/NMF_FDG  
  
##########################################  
#                                        #  
#   Output some useful job information.  #  
#                                        #  
##########################################  
  
# stdout  
echo "------------------------------------------------------"   
echo "SLURM: sbatch is running on ${SLURM_SUBMIT_HOST}"   
echo "SLURM: executing queue is ${SLURM_JOB_PARTITION}"   
echo "SLURM: working directory is ${SLURM_SUBMIT_DIR}"   
echo "SLURM: job identifier is ${SLURM_JOB_ID}"   
echo "SLURM: job name is ${SLURM_JOB_NAME}"   
echo "SLURM: cores per node ${SLURM_CPUS_ON_NODE}"   
echo "SLURM: cores per task ${SLURM_CPUS_PER_TASK}"   
echo "SLURM: node list ${SLURM_JOB_NODELIST}"   
echo "SLURM: number of nodes for job ${SLURM_JOB_NUM_NODES}"   
echo "------------------------------------------------------"  
echo "Command: ${BinDir}/run_extractBasesMT.sh"   
echo "Arguments: Input files: ${inFiles} NumBases: ${numBases} Output directory: ${outputDir}"  
  
( echo -e "Executing in: \c"; pwd )   
( echo -e "Executing at: \c"; date )   
  
# stderr  
echo "------------------------------------------------------" 1>&2  
echo "SLURM: sbatch is running on ${SLURM_SUBMIT_HOST}" 1>&2  
echo "SLURM: executing queue is ${SLURM_JOB_PARTITION}" 1>&2  
echo "SLURM: working directory is ${SLURM_SUBMIT_DIR}" 1>&2  
echo "SLURM: job identifier is ${SLURM_JOB_ID}" 1>&2  
echo "SLURM: job name is ${SLURM_JOB_NAME}" 1>&2  
echo "SLURM: cores per node ${SLURM_CPUS_ON_NODE}" 1>&2   
echo "SLURM: cores per task ${SLURM_CPUS_PER_TASK}" 1>&2  
echo "SLURM: node list ${SLURM_JOB_NODELIST}" 1>&2  
echo "SLURM: number of nodes for job ${SLURM_JOB_NUM_NODES}" 1>&2  
echo "------------------------------------------------------" 1>&2 
echo "Command: ${BinDir}/run_extractBasesMT.sh" 1>&2  
echo "Arguments: Input files: ${inFiles} NumBases: ${numBases} Output directory: ${outputDir}" 1>&2 
  
( echo -e "Executing in: \c"; pwd ) 1>&2  
( echo -e "Executing at: \c"; date ) 1>&2  
  
  
# keep some extra information in the output directory  
mkdir -p ${outputDir}/NumBases${numBases}  
  
# date  
date2save=$(date +"%m-%d-%y")  
  
# output command  
command="${BinDir}/run_extractBasesMT.sh ${MCRROOT} OPNMF ${inFiles} 1 ${numBases} outputDir ${outputDir} saveInterm 1"

echo ${command} > ${outputDir}/NumBases${numBases}/command_${date2save}.txt  
${extractBases} ${MCRROOT} OPNMF ${inFiles} 1 ${numBases} outputDir ${outputDir} saveInterm 1
  
if [ $? != 0 ] ;   
then  
	date_info=$(date)  
	echo "${date_info} : Failure to execute extractBasesMT"  
	echo "${date_info} : Failure to execute extractBasesMT, JobID: ${SLURM_JOB_ID}, method: OPNMF, number of bases:  ${numBases}" >> ${outputDir}/FailedExtractBasesExperimentsMT.txt  

	exit 1  
fi  
