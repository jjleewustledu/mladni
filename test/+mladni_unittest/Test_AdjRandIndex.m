classdef Test_AdjRandIndex < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 09-Jan-2024 22:32:45 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/test/+mladni_unittest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        dataset_home
        testObj
    end
    
    methods (Test)
        function test_ctor(this)
            this.verifyEqual(this.dataset_home, this.testObj.dataset_home)
            disp(this.testObj)
        end
        function test_evaluate_pair_similar(this)
            pwd0 = pushd(this.dataset_home);
            results1 = load(fullfile("NumBases2", "OPNMF", "ResultsExtractBases.mat"));
            results2 = load(fullfile("NumBases4", "OPNMF", "ResultsExtractBases.mat"));
            B1 = results1.B;
            B2 = [results2.B(:,2),results2.B(:,1)];
            %this.view(B1, B2);
            tic
            [ARI,overlap] = this.testObj.evaluate_pair(B1, B2);
            toc
            this.verifyEqual(0.8743428571023, ARI, RelTol=1e-12)
            this.verifyEqual([0.760058717138855;0.934371230466369], overlap, RelTol=1e-12)
            fprintf("%s: ARI = %s\n", stackstr(), mat2str(ARI));
            fprintf("%s: overlap = %s\n", stackstr(), mat2str(overlap));
            popd(pwd0);
        end
        function test_evaluate_pair_identical(this)
            pwd0 = pushd(this.dataset_home);
            results1 = load(fullfile("NumBases2", "OPNMF", "ResultsExtractBases.mat"));
            B1 = results1.B;
            B2 = B1;
            %this.view(B1, B2);
            tic
            [ARI,overlap] = this.testObj.evaluate_pair(B1, B2);
            toc
            this.verifyEqual(1, ARI, RelTol=1e-12)
            this.verifyEqual([1;1], overlap, RelTol=1e-12)
            fprintf("%s: ARI = %s\n", stackstr(), mat2str(ARI));
            fprintf("%s: overlap = %s\n", stackstr(), mat2str(overlap));
            popd(pwd0);
        end
        function test_evaluate_pair_vec(this)
            pwd0 = pushd(this.dataset_home);
            results1 = load(fullfile("NumBases2", "OPNMF", "ResultsExtractBases.mat"));
            B1 = results1.B(:,1);
            B2 = results1.B(:,2);
            %this.view(B1, B2);
            tic
            [ARI,overlap] = this.testObj.evaluate_pair(B1, B2);
            toc
            %this.verifyEqual(1, ARI, RelTol=1e-12)
            %this.verifyEqual([1;1], overlap, RelTol=1e-12)
            fprintf("%s: ARI = %s\n", stackstr(), mat2str(ARI));
            fprintf("%s: overlap = %s\n", stackstr(), mat2str(overlap));
            popd(pwd0);
        end
    end
    
    methods (TestClassSetup)
        function setupAdjRandIndex(this)
            import mladni.*
            this.dataset_home = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "NMF_FDG", "baseline_ad");
            this.testObj_ = AdjRandIndex(this.dataset_home);
        end
    end
    
    methods (TestMethodSetup)
        function setupAdjRandIndexTest(this)
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
        function view(~, B1, B2)
            atl = mlfourd.ImagingContext2(fullfile(getenv("REFDIR"), "MNI152_T1_2mm.nii.gz"));
            ifc1 = copy(atl.imagingFormat);
            ifc1.img = reshape(B1, [91,109,91]);
            ifc2 = copy(atl.imagingFormat);
            ifc2.img = reshape(B2, [91,109,91]);
            atl.view(ifc1, ifc2);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
