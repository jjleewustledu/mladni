classdef NMFRegression3 < handle
    %% line1
    %  line2
    %  
    %  Created 27-Aug-2024 23:03:56 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2689473 (R2024a) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        nmf
        workdir

        cohort_coefficients_sorted
        basis_alluvial
    end

    properties (Dependent)
        N_PATTERNS
        componentsDir
        mask
        niiImgDir
        outputDir
    end
    
    methods  %% GET
        function g = get.N_PATTERNS(this)
            g = this.nmf.N_PATTERNS;
        end
        function g = get.componentsDir(this)
            g = fullfile(this.outputDir, sprintf('NumBases%i', this.N_PATTERNS), 'components');
        end
        function g = get.mask(~)
            g = mlfourd.ImagingContext2( ...
                fullfile(getenv('ADNI_HOME'), 'VolBin', 'mask.nii.gz'));
        end
        function g = get.niiImgDir(this)
            g = fullfile(this.outputDir, sprintf('NumBases%i', this.N_PATTERNS), 'OPNMF', 'niiImg');
        end
        function g = get.outputDir(this)
            g = fullfile(this.workdir, 'baseline_cn');
        end
    end

    methods
        function ic = synthetic_topography(this, opts)
            arguments
                this mladni.NMFRegression3
                opts.param {mustBeText} = "(Intercept)"  % "sexM", "apoe4", "cohortCDR=0,amy+", "cohortCDR=0.5,amy+", "cohortCDR>0.5,amy+"
            end

            T = this.cohort_coefficients_sorted;
            T = T(contains(T.ParamCoefficients, opts.param), :);
            T = natsortrows(T, [], {'ParamCoefficients'});
            encoding = T.Estimate;

            nii = this.basis_alluvial.nifti;
            ifc = copy(nii);
            for b = 1:this.nmf.N_PATTERNS
                ifc.img(:,:,:,b) = encoding(b)*ifc.img(:,:,:,b);
            end
            ifc.img = sum(ifc.img, 4);
            ic = mlfourd.ImagingContext2(ifc);
            ic = ic .* this.mask;
            ic.fileprefix = strrep(strrep(strrep(strrep(opts.param, '(', ''), ')', ''), '>', 'gt'), '=', 'eq');
        end        

        function this = NMFRegression3()
            this.nmf = mladni.NMF();
            this.workdir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');

            ld = load(fullfile(this.workdir, 'CohortCoefficientsSorted.mat'));
            this.cohort_coefficients_sorted = ld.T;
            this.basis_alluvial = mlfourd.ImagingContext2( ...
                fullfile(this.niiImgDir, 'Basis_alluvial.nii'));
        end
    end

    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
