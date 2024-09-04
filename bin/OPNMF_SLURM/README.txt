ToDo:

1. rename submit_dataset_nXX_numBases2to40.sh to something more meaningful that specifies which dataset is used, the number of samples XX, and the range of the number of components that were evaluated (here 2 to 40 with step 2)

2. go in submit_dataset_nXX_numBases2to40.sh and change the following:
	a) inFiles="path to .csv file" -- this variable should contain the full path to a .csv file that contains in every row the full path to each nifti file 
	b) outputDir="OutputDirectory" -- replace OutputDirectory with the path to the directory
	c) scriptDir="path to VolBin" -- change the variabe to the path that points to the VolBin that I shared with you

3. "for b in 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40" -- change the values here to the include the number(s) of the components you want to estimate. Based on our discussion, it might make sense to run with a step of 1 and a max of 20. In that case, replace line with "for b in {1..20}"

4. in VolBin, open the create_submit_extractBasesMT.sh and change the following:
	a) echo '#PBS -M wustlid@wustl.edu  ' -- replace with your email
	b) echo '### Send email only if the job is aborted (flag a), when the job starts (flag b), or ends (flag e)  '
	   echo '#PBS -m a
	
	Currently, one get an email only when the job is aborted. You can change the flags to get an email when the job starts or ends.

	c) Change the requested resources if necessary 
	
	echo '### Select 1 nodes with 4 CPU  '
	echo '#PBS -l nodes=1:ppn=4  '
	echo '### Select 24Gb RAM  '
	echo '#PBS -l mem=24GB  '
	echo '### Select 1h wall time  '
	echo '#PBS -l walltime=24:00:00  '

	Currently, the script asks for 24Gb RAM, 4 CPUs and 24h of compute time. Asking for more CPUs will make the executable finish faster. The more resources one asks, the more time it will probably take for the job to start. If an executable needs more resources than the ones requested, it will fail. I am not sure about the exact requirements here as this depends both on the number of samples and the dimensionality of data. I am almost certain that the nmf computations will take more than 24h (especially for the higher number of components). However, I would be sKeptical about requesting more time as it seems that the queueing system heavily penalizes longer jobs (small number of max running jobs). Note that the code saves intermediate results and can restart from these intermediate points by simply running the same command and pointing to the same output directory.

5) in: echo 'BinDir="path to VolBin/extractBasesMT"  ' change the path to CortThickBin with the path of the directory where you are going to move the code I shared with you

6) run the renamed submit_dataset_nXX_numBases2to40.sh to submit the nmf experiments to the job scheduler. This is done by typing the command ./submit_dataset_nXX_numBases2to40.sh
