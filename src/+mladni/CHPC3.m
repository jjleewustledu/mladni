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
        function propcluster_tiny()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '10000'; % in MB; deepmrseg requires 10 GB; else 6 GB
            c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = 'test';
            c.AdditionalProperties.WallTime = '4:00:00'; % 4 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function propcluster()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '10000'; % in MB; deepmrseg requires 10 GB; else 6 GB
            c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = '';
            c.AdditionalProperties.WallTime = '24:00:00'; % 24 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function setenvs()
            [~,r] = system('hostname');
            if ~contains(r, 'cluster')
                return
            end

            setenv('TMPDIR', '/scratch/jjlee/tmp') % worker nodesk

            setenv('SINGULARITY_HOME', '/scratch/jjlee/Singularity')
            setenv('ADNI_HOME', '/scratch/jjlee/Singularity/ADNI') 
            setenv('AFNIPATH', '/export/afni/afni-20.3.03/linux_openmp_64')
            setenv('ANTSPATH', '/export/ants/ants-2.3.5/bin')
            setenv('DEBUG', '');
            setenv('FREESURFER_HOME', '/home/jjlee/.local/freesurfer/freesurfer-7.3.2')
            setenv('FSLDIR', '/export/fsl/fsl-6.0.5')

            setenv('FSLOUTPUTTYPE', 'NIFTI_GZ')
            setenv('FSLMULTIFILEQUIT', 'TRUE')
            setenv('FSLTCLSH', fullfile(getenv('FSLDIR'),'bin','fsltclsh'))
            setenv('FSLWISH', fullfile(getenv('FSLDIR'),'bin','fslwish'))
            setenv('FSLLOCKDIR', '')
            setenv('FSLMACHINELIST', '')
            setenv('FSLREMOTECALL', '')
            setenv('FSLREMOTECALL', 'cuda.q')
            setenv('PYOPENGL_PLATFORM', 'osmesa')

            setenv('REFDIR', '/home/jjlee/.local/atlas')
            setenv('RELEASE', '/home/jjlee/.local/lin64-tools')            
            setenv('PATH', ...
                strcat(getenv('RELEASE'), ':', ...
                       getenv('AFNIPATH'), ':', ...
                       fullfile(getenv('FREESURFER_HOME'), 'bin'), ':', ...
                       fullfile(getenv('FSLDIR'), 'bin'), ':', ...
                       '/export/singularity/singularity-3.9.0/bin', ':', ...
                       getenv('PATH')))
            setenv('LD_LIBRARY_PATH', ...
                strcat('/usr/lib64', ':', getenv('LD_LIBRARY_PATH'))) % need libOSMesa.so.8 for fsleyes render
                   
            %disp("mladni.CHPC3.setenvs():getenv('PATH'):")
            %disp(getenv('PATH'))
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
