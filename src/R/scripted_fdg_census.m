% ensure ~jjlee/Documents/MATLAB/MatlabRegistry.m in Matlab path
registry = MatlabRegistry.instance();

% set environment
%setenv('SINGULARITY_HOME', '/scratch/jjlee/Singularity');
%setenv('ADNI_HOME', fullfile(getenv(SINGULARITY_HOME), 'ADNI'));

ad = mladni.AdniDemographics(); % ~jjlee/MATLAB-Drive/mladni/scr/+mladni/AdniDemographics

% uncomment to see truncated master table of all available FDG data from ADNI/LONI through 2021-nov-23.
%disp(ad.table_fdg1) 

% remove trailing semi-colons to examine truncated tables in command window
ad.table_cn;
ad.table_preclinical;
ad.table_cdr_0p5_apos_emci;
ad.table_cdr_0p5_apos_lmci;
ad.table_cdr_gt_0p5_apos;

