classdef CHPC3
    %% Provides support functions for using the Matlab Parallel Server at CHPC3.
    %  
    %  Created 07-Apr-2022 16:13:53 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function clean_tempdir()
            try
                deleteExisting(fullfile(tempdir, '*.nii*'));
                deleteExisting(fullfile(tempdir, '*.4dfp.*'));
            catch ME
                disp(ME)
            end
        end
        function setenvs()
            setenv('TMPDIR', '/scratch/jjlee/tmp') % worker nodes

            setenv('ADNI_HOME', '/home/aris_data/ADNI_FDG') 
            setenv('AFNIPATH', '/export/afni/afni-20.3.03/linux_openmp_64')
            setenv('ANTSPATH', '/export/ants/ants-2.3.5/bin')
            setenv('DEBUG', '');
            setenv('FREESURFER_HOME', '/export/freesurfer/freesurfer-7.2.0')
            setenv('FSLDIR', '/export/fsl/fsl-6.0.5')

            setenv('FSLOUTPUTTYPE', 'NIFTI_GZ')
            setenv('FSLMULTIFILEQUIT', 'TRUE')
            setenv('FSLMULTIFILEQUIT', 'TRUE')
            setenv('FSLTCLSH', fullfile(getenv('FSLDIR'),'bin','fsltclsh'))
            setenv('FSLWISH', fullfile(getenv('FSLDIR'),'bin','fslwish'))
            setenv('FSLLOCKDIR', '')
            setenv('FSLMACHINELIST', '')
            setenv('FSLREMOTECALL', '')
            setenv('FSLREMOTECALL', 'cuda.q')

            setenv('REFDIR', '/home/aris_data/ADNI_FDG/atlas')
            setenv('RELEASE', '/home/aris_data/ADNI_FDG/lin64-tools')            
            setenv('PATH', ...
                strcat(getenv('RELEASE'), ':', ...
                       getenv('AFNIPATH'), ':', ...
                       fullfile(getenv('FREESURFER_HOME'), 'bin'), ':', ...
                       fullfile(getenv('FSLDIR'), 'bin'), ':', ...
                       '/export/gsl/gsl-2.7.1/bin:/export/afni/afni-20.3.03/linux_openmp_64:/export/singularity/singularity-3.9.0/bin:/export/cuda/cuda-10.2/bin:/usr/bin:/export/freesurfer/freesurfer-7.2.0/bin:/export/freesurfer/freesurfer-7.2.0/fsfast/bin:/export/freesurfer/freesurfer-7.2.0/tktools:/export/freesurfer/freesurfer-7.2.0/mni/bin:/export/ants/ants-2.3.5/bin:/export/fsl/fsl-6.0.5/bin:/home/aris_data/ADNI_FDG/lin64-tools:/home/jjlee/.local/bin:/home/jjlee/bin:/usr/local/bin:/usr/local/sbin:/usr/sbin'), ':', ...
                       fullfile('/export/singularity/singularity-3.9.0/bin'))
                   
            %disp("mladni.CHPC3.setenvs():getenv('PATH'):")
            %disp(getenv('PATH'))
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
