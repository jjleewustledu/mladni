classdef Test_NMF < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 08-Jul-2023 13:03:54 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/test/+mladni_unittest.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        N_PATTERNS = mladni.NMF.N_PATTERNS
    end

    properties
        study_design = "cross-sectional"
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mladni.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_prep_regressions(this)
        end
        function test_table_covariates(this)
            % mladni.AdniDemographics.table_covariates must be well-formed with pattern-averaged FDG SUVR metrics

            nmfcov_long = mladni.NMFCovariates(study_design="longitudinal");
            Tc_long = nmfcov_long.table_covariates %#ok<NOPRT> % visual inspect
            figure; plot(Tc_long.Components)
            this.verifyEqual(size(Tc_long), [3415, 135+4]) % Dlicv, PVE1, RegErr, Components
            this.verifyEqual(size(Tc_long.Components), [3415, this.N_PATTERNS])

            nmfcov_cs = mladni.NMFCovariates(study_design="cross-sectional");
            Tc_cs = nmfcov_cs.table_covariates %#ok<NOPRT> % visual inspect
            figure; plot(Tc_cs.Components)
            this.verifyEqual(size(Tc_cs), [3415, 135+4]) % Dlicv, PVE1, RegErr, Components
            this.verifyEqual(size(Tc_cs.Components), [3415, this.N_PATTERNS])
        end
        function test_table_dlicv(this)
            nmfcov = mladni.NMFCovariates(study_design="longitudinal");
            t = nmfcov.table_dlicv % visual inspect
            plot(t.Dlicv)
        end
        function test_table_fdg4(this)
            % mladni.AdniDemographics.table_fdg4 must be well-formed with fqfn for registered FDG

            ad = mladni.AdniDemographics(study_design="longitudinal");
            T3 = ad.table_fdg3;
            this.verifyEqual(size(T3), [3478, 134])
            T4 = ad.table_fdg4;
            this.verifyEqual(size(T4), [3478, 135])

            ad = mladni.AdniDemographics(study_design="cross-sectional");
            T3 = ad.table_fdg3;
            this.verifyEqual(size(T3), [3478, 134])
            T4 = ad.table_fdg4;
            this.verifyEqual(size(T4), [3478, 135])
        end
        function test_table_filelist(this)
            nmfcov = mladni.NMFCovariates(study_design="longitudinal");
            nmfcov.table_filelist % visual inspect
        end
        function test_table_imagedataIDs(this)
            nmfcov = mladni.NMFCovariates(study_design="longitudinal");
            nmfcov.table_imagedataIDs % visual inspect
        end
        function test_table_pve1(this)
            nmfcov = mladni.NMFCovariates(study_design="longitudinal");
            t = nmfcov.table_pve1 % visual inspect
            plot(t.PVE1)
        end
        function test_table_regerr(this)
            nmfcov = mladni.NMFCovariates(study_design="longitudinal");
            t = nmfcov.table_regerr % visual inspect
            N_nan = sum(isnan(t.RegErr));
            fprintf("N_nan -> %g\n", N_nan)
            figure; histogram(t.RegErr)
            figure; plot(t.RegErr)
        end
        function test_table_selectedComponentWeightedAverageNIFTI(this)
            %% Compute expensive.
            %  table_fdg4 ~ 3478 rows; table_selectedComponentWeightedAverageNIFTI ~ 3470 rows,
            %  because of misregistrations that prevented construction of *_pet_on_T1w_Warped_dlicv.nii.gz.

            nmfcov_cs = mladni.NMFCovariates(study_design="cross-sectional");
            t_cs = nmfcov_cs.table_selectedComponentWeightedAverageNIFTI % visually inspect
            figure; plot(t_cs.Components)

            nmfcov_long = mladni.NMFCovariates(study_design="longitudinal");
            t_long = nmfcov_long.table_selectedComponentWeightedAverageNIFTI % visually inspect
            figure; plot(t_long.Components)

            this.verifyEqual(all(size(t_cs), size(t_long)))
        end
    end
    
    methods (TestClassSetup)
        function setupNMF(this)
            import mladni.*
            this.testObj_ = NMF(study_design=this.study_design);
        end
    end
    
    methods (TestMethodSetup)
        function setupNMFTest(this)
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
