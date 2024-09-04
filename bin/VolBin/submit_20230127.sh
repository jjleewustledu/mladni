#!/bin/bash
launch=/scratch/jjlee/Singularity/ADNI/VolBin
${launch}/submit_nmf_dataset.sh -d baseline_cn_repA -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cn_repB -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_preclinical_repA -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_preclinical_repB -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_gt_0p5_apos_repA -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_gt_0p5_apos_repB -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_repA -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_apos_repB -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_aneg_repA -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_aneg_repB -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_anan_repA -k ${launch}/mask.nii.gz
${launch}/submit_nmf_dataset.sh -d baseline_cdr_0p5_anan_repB -k ${launch}/mask.nii.gz
