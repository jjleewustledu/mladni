classdef Test_AdniMerge < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 02-Mar-2023 11:55:03 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/test/+mladni_unittest.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mladni.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_ad(this)
            ad = mladni.AdniDemographics;
            ad.table_fdg
        end
        function test_table_fdg(this)
            ad = mladni.AdniDemographics;
            ad.table_fdg
        end
        function test_table_fdg1(this)
            % no memory cache
            ad = mladni.AdniDemographics;
            % no filesystem cache
            deleteExisting( ...
                fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', 'AdniDemographics_table_fdg1.mat'))

            % build cache
            ad.table_fdg1
        end
        function test_am(this)
            am = mladni.AdniMerge;
            am.table_merge
        end
    end
    
    methods (TestClassSetup)
        function setupAdniMerge(this)
            import mladni.*
            this.testObj_ = AdniMerge();
        end
    end
    
    methods (TestMethodSetup)
        function setupAdniMergeTest(this)
            this.testObj = this.testObj_;
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
