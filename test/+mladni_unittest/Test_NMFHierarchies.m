classdef Test_NMFHierarchies < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 20-Feb-2024 12:55:36 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/test/+mladni_unittest.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
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
        function test_ctor(this)
            this.verifyTrue(isfolder(this.testObj.home))
        end
        function test_table_patt_weighted_fdg(this)
            T = this.testObj.table_patt_weighted_fdg(N_bases_target=24)
        end
    end
    
    methods (TestClassSetup)
        function setupNMFHierarchies(this)
            import mladni.*
            this.testObj_ = NMFHierarchies();
        end
    end
    
    methods (TestMethodSetup)
        function setupNMFHierarchiesTest(this)
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
