# NMF SLURM

This repository contains code for running Non-negative Matrix Factorization (NMF) on the CHPC3 with SLURM.  It contains compiled MATLAB code for NMF ([`nmf_slurm/VolBin/extractBasesMT`](nmf_slurm/VolBin/extractBasesMT)).  This readme will describe the edits needed to run the code for new data.  Most instructions are copied directly from the previously existing README.txt (also available in this repo).

## Instructions

1. Rename `submit_dataset_nXX_numBases2to40.sh` to something more meaningful that specifies which dataset is used, the number of samples XX, and the range of the number of components that were evaluated (here 2 to 40 with step 2).

2. Go in `submit_dataset_nXX_numBases2to40.sh` and change the following:

   - `inFiles="/full/path/to/csv/file"` -- this variable should contain the full path to a .csv file that contains in every row the full path to each nifti file 
   - `outputDir="/full/path/to/output_directory"` -- replace OutputDirectory with the full path to the directory
   - `scriptDir="/full/path/to/VolBin"` -- change the variabe to the full path that points to the VolBin folder
   - `ranks="2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40"` -- change the values here to the include the number(s) of the components you want to estimate.

3. In VolBin, open the`create_submit_extractBasesMT.sh` and change the following:

   - `PARAMETERS=(saveInterm 1)` -- add any additional parameters for NMF here.  For a list of valid parameters, see the included `parseInput_REFERENCE.m` file.  Note that in addition to the parameters specified there, **you can add a mask here** by passing  `mask path/to/mask.nii.gz`.  This should be a mask image with same dimensions as your input images.  

   - Set email preferences.  *Note that by default I have commented out emailing, but this can be added back in*.

     - `echo '#SLURM --mail-user=wustlid@wustl.edu  ' ` -- replace with your email.  
     - `echo '#SBATCH --mail-type=END,FAIL,TIME_LIMIT'` -- set the scenarios which trigger an email.  [See documentation for sbatch](https://slurm.schedmd.com/sbatch.html#OPT_mail-type).

   - Change the requested resources if necessary.  There is a note included below which discusses how to pick resource requests.  Furthermore, [see documentation for sbatch](https://slurm.schedmd.com/sbatch.html#OPT_mail-type):

     ```
     echo '### Start job with priority - default is zero  '
     echo '#SLURM --priority=0  '
     echo '### Select 1 nodes with 4 CPU  '
     echo '#SLURM --nodes=1  '
     echo '#SLURM --ntasks-per-node=4  '
     echo '### Select 24Gb RAM  '
     echo '#SLURM --mem=128G  '
     echo '### Select wall time  '
     echo '#SLURM -t 24:00:00  '
     ```

4. Run the renamed `submit_dataset_nXX_numBases2to40.sh` to submit the NMF experiments to the job scheduler. This is done by typing the command: `./submit_dataset_nXX_numBases2to40.sh`. 

## Resources

This note was included in the original README, discussing resources:

> Currently, the script asks for 24Gb RAM, 4 CPUs and 24h of compute time. Asking for more CPUs will make the executable finish faster. The more resources one asks, the more time it will probably take for the job to start. If an executable needs more resources than the ones requested, it will fail. I am not sure about the exact requirements here as this depends both on the number of samples and the dimensionality of data. I am almost certain that the nmf computations will take more than 24h (especially for the higher number of components). However, I would be sKeptical about requesting more time as it seems that the queueing system heavily penalizes longer jobs (small number of max running jobs). Note that the code saves intermediate results and can restart from these intermediate points by simply running the same command and pointing to the same output directory.

## Analyses

There is also some of Aris' MATLAB code for NMF analysis, which is included in the `extra_scripts` folder.  These can be used for analyzing the NMF results.  Some highlights:

- `calcRecError.m`: plots curves of the reconstruction error (difference between input data matrix and the product of the factorized outputs) for different ranks.
- `evaluateReproducibility.m`: compares the similarity of components output by two different NMF runs with the same ranks.  Typical use of this would involve splitting your sample into two matched subsamples.
- `calculateComponentWeightedAverageNIFTI.m`: projects images onto the learned components for a given rank, and provides the weighted loading of each component for the input.  This can be used to create loadings with meaningful units.

Note, these scripts depend on [this MATLAB package](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), which is included in the `NIfTI_20140122` folder.