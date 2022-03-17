classdef Test_FDG < matlab.unittest.TestCase
    %% These unit tests for class mladni.FDG serve purposes of documentation as well as testing.
    %  While confined for use within Matlab, mature code bases may be deployed as Matlab compiled objects or as Docker 
    %  containers.  
    %
    %  E.g.:
    %  >> runtests("/Users/jjlee/MATLAB-Drive/mladni/test/+mladni_unittest/Test_FDG.m","ProcedureName","test_flirt")
    %  
    %  Created 11-Mar-2022 21:02:05 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/test/+mladni_unittest.
    %  Developed on Matlab 9.11.0.1873467 (R2021b) Update 3 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function test_sercluster()
            % Get a handle to the cluster
            c = parcluster
            c.AdditionalProperties

            % Submit job to query where MATLAB is running on the cluster
            j = c.batch(@pwd, 1, {}, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false)

            % Query job for state
            for ip = 1:10
                pause(60)
                if strcmp(j.State, 'finished') || strcmp(j.State, 'failure')
                    break
                end
            end

            try
                % If state is finished, fetch the results
                j.fetchOutputs{:}
                % Delete the job after results are no longer needed
                %j.delete
            catch ME
                handwarning(ME)
            end            
            j.getTaskSchedulerIDs()
            c.getDebugLog(j.Tasks(end))
        end        
        function test_parcluster()
            % Get a handle to the cluster
            c = parcluster;
            c.AdditionalProperties

            % Submit a batch pool job using 4 workers for 16 simulations
            j = c.batch(@mladni_unittest.Test_FDG.parallel_example, 1, {16}, 'Pool', 4, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);

            % View current job status
            for ip = 1:10
                pause(60)
                if strcmp(j.State, 'finished') || strcmp(j.State, 'failure')
                    break
                end
            end
            
            try
                % Fetch the results after a finished state is retrieved
                j.fetchOutputs{:}
                % Delete the job after results are no longer needed
                %j.delete
            catch ME
                handwarning(ME)
            end   
            j.getTaskSchedulerIDs()
            c.getDebugLog(j)
            disp('see also /mnt/beegfs/home/jjlee/RESULTS.mat')
        end        
        function test_sercluster_foo()
            % Get a handle to the cluster
            c = parcluster
            c.AdditionalProperties

            % Submit job to query where MATLAB is running on the cluster
            j = c.batch(@foo, 1, {}, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false)

            % Query job for state
            for ip = 1:10
                pause(60)
                if strcmp(j.State, 'finished') || strcmp(j.State, 'failure')
                    break
                end
            end

            try
                % If state is finished, fetch the results
                j.fetchOutputs{:}
                % Delete the job after results are no longer needed
                %j.delete
            catch ME
                handwarning(ME)
            end            
            j.getTaskSchedulerIDs()
            c.getDebugLog(j.Tasks(end))
            
            function g = foo()    
                
                %mlbash('module restore')
                setenv('ADNI_HOME', '/home/aris_data/ADNI_FDG')
                setenv('ANTSPATH', '/export/ants/ants-2.3.5/bin')
                setenv('DEBUG', '');
                setenv('FREESURFER_HOME', '/export/freesurfer/freesurfer-7.2.0')
                setenv('FSLDIR', '/export/fsl/fsl-6.0.5')
                setenv('RELEASE', '/home/aris_data/ADNI_FDG/lin64-tools')
                
                setenv('PATH', ...
                    strcat(getenv('RELEASE'), ':', ...
                           fullfile(getenv('FREESURFER_HOME'), 'bin'), ':', ...
                           fullfile(getenv('FSLDIR'), 'bin'), ':', ...
                           getenv('PATH')))            

                [~,g] = system('fslorient');
            end
        end
        function test_parcluster_foo()
            % Get a handle to the cluster
            c = parcluster;
            c.AdditionalProperties.WallTime = '02:00:00';
            c.AdditionalProperties.Node = 4;
            c.AdditionalProperties.MemUsage = '10000';
            disp(c.AdditionalProperties)

            % Submit a batch pool job using 4 workers for 16 simulations
            j = c.batch(@foo, 1, {4}, 'Pool', 2, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);

            % View current job status
            for ip = 1:30
                pause(60)
                if strcmp(j.State, 'finished') || strcmp(j.State, 'failure')
                    break
                end
            end
            
            try
                % Fetch the results after a finished state is retrieved
                j.fetchOutputs{:}
                % Delete the job after results are no longer needed
                %j.delete
            catch ME
                handwarning(ME)
            end   
            j.getTaskSchedulerIDs()
            c.getDebugLog(j)
            disp('see also /mnt/beegfs/home/jjlee/')            
            
            function g = foo(iter)    
                g = cell(1,iter);
                parfor idx = 1:iter

                    %mlbash('module restore')
                    setenv('ADNI_HOME', '/home/aris_data/ADNI_FDG')
                    setenv('ANTSPATH', '/export/ants/ants-2.3.5/bin')
                    setenv('DEBUG', '');
                    setenv('FREESURFER_HOME', '/export/freesurfer/freesurfer-7.2.0')
                    setenv('FSLDIR', '/export/fsl/fsl-6.0.5')
                    setenv('RELEASE', '/home/aris_data/ADNI_FDG/lin64-tools')

                    setenv('PATH', ...
                        strcat(getenv('RELEASE'), ':', ...
                               fullfile(getenv('FREESURFER_HOME'), 'bin'), ':', ...
                               fullfile(getenv('FSLDIR'), 'bin'), ':', ...
                               getenv('PATH')))  

                    [~,g{idx}] = system('fslorient');
                end
            end
        end 
        function [t,A] = parallel_example(iter)
            if nargin==0
                iter = 8;
            end
            
            disp('Start sim')
            
            t0 = tic;
            parfor idx = 1:iter
                A(idx) = idx;
                pause(2)
                idx
            end
            t = toc(t0);
            
            disp('Sim completed')
            
            save('RESULTS.mat', 'A')            
        end
    end
    
    properties
        fdg
        proc = 'CASU'
        t1w
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mladni.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_ctor(this)
            disp(this.testObj)
        end
        function test_call(this)
            this.testObj.call();
            this.testObj.view();
            this.testObj.finalize();
        end
        function test_call_site941(this)
            setenv('DEBUG', '1')
            globbed = globT( ...
                fullfile(getenv('ADNI_HOME'), sprintf('bids/rawdata/sub-941S*/ses-*/pet/*-%s_pet.nii.gz', this.proc)));
            %len = length(globbed);
            %globbed = globbed(rand(1,len) > 0.9);
            fprintf('test_call_site941:  sampling %i FDGs\n', length(globbed)); % ~60 @ 185 sec
            for ig = 1:length(globbed) % parfor fails for bids on pascal:/scratch
                fdg_ = mlfourd.ImagingContext2(globbed{ig});
                obj = mladni.FDG(fdg_);
                obj.call();
                obj.finalize();
            end
        end
        function test_call_resolve_site941(this)
            setenv('DEBUG', '1')
            globbed = globT( ...
                fullfile(getenv('ADNI_HOME'), sprintf('bids/rawdata/sub-941S*/ses-*/pet/*-%s_pet.nii.gz', this.proc)));
            %len = length(globbed);
            %globbed = globbed(rand(1,len) > 0.9);
            fprintf('test_call_site941:  sampling %i FDGs\n', length(globbed)); % ~60 @ 185 sec
            for ig = 1:length(globbed) % parfor fails for bids on pascal:/scratch
                %pwd0 = pushd(myfileparts(globbed{ig}));
                fdg_ = mlfourd.ImagingContext2(globbed{ig});
                obj = mladni.FDG(fdg_);
                obj.call_resolve();
                obj.finalize();
                %popd(pwd0);
            end
        end
    end
    
    methods (TestClassSetup)
        function setupFDG(this)
            import mladni.*

            setenv('ADNI_HOME', fullfile(getenv('SINGULARITY_HOME'), 'ADNI', ''));
            this.fdg = mlfourd.ImagingContext2( ...
                fullfile(getenv('ADNI_HOME'), ...
                'bids/rawdata/sub-022S0543/ses-20060619/pet', ...
                strcat('sub-022S0543_ses-20060619100038_trc-FDG_proc-', this.proc, '_pet.nii.gz')));
            this.testObj_ = FDG(this.fdg);
        end
    end
    
    methods (TestMethodSetup)
        function setupFDGTest(this)
            if ishandle(this.testObj)
                this.testObj = copy(this.testObj_);
            else
                this.testObj = this.testObj_;
            end
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
