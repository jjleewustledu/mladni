# mladni
mladni supports machine learning models for neuroimaging, emphasizing studies of Alzheimer's Dementia and using data from [the Alzheimer's Disease Neuroimaging Initiative](http://adni.loni.usc.edu/about/).  

Existing and planned software dependencies include [Matlab](https://www.mathworks.com/), [TensorFlow](https://www.tensorflow.org), and [MONAI](https://monai.io/).  Some software dependencies are described in [public repositories by John Lee](https://github.com/jjleewustledu?tab=repositories).

# Getting Started

## Installation

This work is in progress.  Detailed requirements are enumerated in [mladni_unittest.requiredFilesAndProducts.m](https://github.com/sotiraslab/mladni/blob/main/test/%2Bmladni_unittest/requiredFilesAndProducts.m), available from [public repositories by John Lee](https://github.com/jjleewustledu?tab=repositories) and [Matlab](https://www.mathworks.com/).

## Use case:  explore the database of previously curated ADNI data.

ADNI comprises a consortium of research entities which have collected clinical and imaging data from human populations at risk of developing Alzheimer's Disease.   [ADNIMERGE](https://adni.bitbucket.io/) is the entity with the most comprehensive curation of data, but it is significantly constrained by legacy conventions since the beginnings of the consortium in 2004.  This project, mladni, avoids the R-based functionality of ADNIMERGE, preferring instead to use files of flat comma-separated values (csv).  These files of csv may be structured for new research aims, retaining only weak dependencies upon legacy conventions of ADNIMERGE.  See mladni.AdniMerge and mladni.AdniDemographics for examples of accessing ADNIMERGE with weak dependencies.  mladni.AdniDemographics.table_fdg1() is a working example of managing heterogeneous ADNIMERGE sources for a project focused on PET with FDG. 

## Use case:  create bids filesystem hierarchies for neuroimaging downloaded from ADNI.

mladni.{[AdniBids](https://github.com/sotiraslab/mladni/blob/main/src/%2Bmladni/AdniBids.m),[AdniBidsT1w](https://github.com/sotiraslab/mladni/blob/main/src/%2Bmladni/AdniBidsT1w.m),[AdniBidsFdg](https://github.com/sotiraslab/mladni/blob/main/src/%2Bmladni/AdniBidsFdg.m)} provide working examples.  See static methods create() for ideas that have been useful to date.

## Use case:  prepare T1-weighted imaging and PET FDG for training a machine learning model such as non-negative matrix factorization.  

Models which depend on linear mappings between training data and predictive parameters often require careful co-registration of neuroimaging within a standardized atlas space.  A common workflow is:

1. Co-register PET to T1w anatomy.
   - ADNI provides PET with various [types of processing](http://adni.loni.usc.edu/methods/pet-analysis-method/pet-analysis/).  Processing that includes co-registration of dynamic frames, averaging frames, standardizing image and voxel size, and imposing uniform resolution, or CASU, will be useful for most projects.
   - PET and T1w anatomy have distinct features.  The primary correspondence between these modalities is anatomic geometry.  For this reason, use rigid body registration.
   - The anatomical geometry of structures inside the skull are largely unchanging for short periods of time.
     + Co-register PET with T1w anatomy nearest in time.
     + Acute medical or surgical events can produce gross alterations of anatomy correspondence.
   - The anatomical geometry of the neck can vary significantly between scans, especially between scanners.
     + Use brain masks as weights for T1w.  See documentation from [FSL FLIRT](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/UserGuide).
   - PET scans have structured scattering phenomena outside of the skull that can disrupt co-registration.  
     + Use generous masks with thresholds ~ 0.0001*max(PET), or no masks, as weights for PET.  
   - PET scanners produce NIfTI objects which are most always in radiological and left-anterior-superior orientation.  Many MRI scanners produce NIfTI in neurological and right-anterior-superior orientation.
     + Prefer using [AFNI's](https://afni.nimh.nih.gov/afni/doc/help/3dresample.html) "3dresample".  [FSL's](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Orientation%20Explained) "fslorient", "fslswapdim", and "reorient2std" are much more difficult to use, and should not be attempted without extensive testing.
   - See [src/mladni.FDG](https://github.com/sotiraslab/mladni/blob/main/src/%2Bmladni/FDG.m).{flirt_t1w2atl,flirt_fdg2t1w2atl} for a worked example.
   - See [doc/"2022 Mar 11.pdf"](https://github.com/sotiraslab/mladni/blob/main/doc/2022%20Mar%2011.pdf) for more information.

2. Warp PET and T1w anatomy to a standard atlas.
   - See [ANTs](https://github.com/ANTsX/ANTs/wiki/Forward-and-inverse-warps-for-warping-images,-pointsets-and-Jacobians).
   - See [mladni.FDG](https://github.com/sotiraslab/mladni/blob/main/src/%2Bmladni/FDG.m).CreateJacobianDeterminantImage.

3. Retain only neuroimaging from the brain, cerebellum and brainstem.
   - Consider using FSL's standardized brain masks, such as $FSL/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz.
   - For some use cases, [DeepMRSeg](https://github.com/CBICA/DeepMRSeg) may be appropriate.

4. Ensure robust processing for large quantities of neuroimaging, e.g., ~4000 PET FDG for ADNI.
   - Carefully visualize intermediates while developing the workflow, e.g., for ~10 imaging sessions.
   - Visualize simplified outputs of workflows, such as co-registered, warped neuroimaging, e.g., for ~100 imaging sessions.
   - Produce quality control metrics for assessing all neuroimaging.
     + Visualize the distribution of final cost-functions following co-registrations, e.g., normalized mutual information.
     + Generate maps showing the average of co-registered neuroimaging.
     + Generate maps showing the variance of co-registered neuroimaging.
     + See also similar functionality from [mcflirt](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MCFLIRT).

For a worked example, see [mladni.FDG.call()](https://github.com/sotiraslab/mladni/blob/main/src/%2Bmladni/FDG.m#CALL).

A working implementations of Matlab distributed computing is [mladni.FDG.parcluster()](https://github.com/sotiraslab/mladni/blob/main/src/%2Bmladni/FDG.m#PARCLUSTER).
