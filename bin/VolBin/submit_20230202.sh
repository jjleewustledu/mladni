#!/bin/bash
launch=/scratch/jjlee/Singularity/ADNI/VolBin

#${launch}/submit_nmf_dataset.sh -d baseline_cn -k ${launch}/mask.nii.gz
#exit





${launch}/submit_nmf_dataset.sh -d baseline_preclinical -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_gt_0p5_apos -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_emci -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_mci -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_lmci -k ${launch}/mask.nii.gz

for i in {1..20}
do
${launch}/submit_nmf_dataset.sh -d baseline_cn_repA${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cn_repB${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_preclinical_repA${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_preclinical_repB${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_gt_0p5_apos_repA${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_gt_0p5_apos_repB${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_emci_repA${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_emci_repB${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_mci_repA${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_mci_repB${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_lmci_repA${i} -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_lmci_repB${i} -k ${launch}/mask.nii.gz
done



