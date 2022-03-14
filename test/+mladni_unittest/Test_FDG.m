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
    
    properties
        fdg
        proc = 'CAS'
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
