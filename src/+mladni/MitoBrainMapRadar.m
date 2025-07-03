classdef MitoBrainMapRadar < handle & mladni.NMFRadar
    %% Use 24-rank PoMC to generate radar plots of continuous pattern samples from ADNI/MitochondrialProfile.
    %  
    %  Created 02-Apr-2025 13:29:57 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.2.0.2863752 (R2024b) Update 5 for MACA64.  Copyright 2025 John J. Lee.
    
    properties
        intracranial_mask
        mito_fps
        mito_profile_home
        mito_profiles
        basis_all
        basis_max
    end

    methods

        function pr = pattern2prob(this, ifc, p)
            ifc_ = copy(ifc);
            ifc_.img = ifc.img(:,:,:,p);
            ic = mlfourd.ImagingContext2(ifc_);
            ic = ic .* this.intracranial_mask;
            ic = ic / this.basis_max;
            pr = ic.imagingFormat.img;
        end

        function T = table(this, varargin)
            %% Table T has variables:
            %  ProfileLabel, Estimate, LowerCI, UpperCI for MitoBrainMap

            if ~isempty(this.table_mito_brain_map_radar_)
                T = this.table_mito_brain_map_radar_;
                return
            end

            table_fqfn = fullfile(this.mito_profile_home, stackstr() + ".mat");
            if isfile(table_fqfn)
                ld = load(table_fqfn);
                T = ld.T;
                return
            end

            ProfileLabel = [];
            Estimate = [];
            LowerCI = [];
            UpperCI = [];
            patterns_ifc = mlfourd.ImagingFormatContext2(this.basis_all);

            for mfp = this.mito_fps
                mfp_ifc = mlfourd.ImagingFormatContext2( ...
                    fullfile(this.mito_profile_home, mfp + ".nii.gz") ...
                );
                tic
                fprintf("\n%s: %s started\n", stackstr(), mfp);
                fprintf("\tp -> ");
                for p = 1:this.N_patterns
                    fprintf("%i ", p);
                    ProfileLabel = [ProfileLabel; mfp]; %#ok<*AGROW>

                    prob = this.pattern2prob(patterns_ifc, p);
                    [e, lci, uci] = this.bootstrapROIMean( ...
                        mfp_ifc.img, ...
                        prob ...
                    );
                    Estimate = [Estimate; e];
                    LowerCI = [LowerCI; lci];
                    UpperCI = [UpperCI; uci];
                end
                toc
            end

            T = table(ProfileLabel, Estimate, LowerCI, UpperCI);
            save(table_fqfn, "T");
            this.table_mito_brain_map_radar_ = T;
        end

        function T = table_intercept(this, varargin)
            matfile = 'CohortCoefficientsSorted.mat';
            ld = load(fullfile(this.num_bases_dir, matfile));
            PC = ld.T;
            intercept = PC(contains(PC.ParamCoefficients, "Intercept"), :);

            Estimate = ascol(intercept.Estimate);
            LowerCI = 1.96 * ascol(intercept.StdErr);
            UpperCI = LowerCI;

            T = table(Estimate, LowerCI, UpperCI);
        end

        function T = table_over_intercept(this, varargin)
            %% Table T has variables:
            %  ProfileLabel, Estimate, LowerCI, UpperCI for MitoBrainMap,
            %  but the variables describe those from table() normalized by the SUVR
            %  intercept provided by mladni.MitoBrainMapRadar.table_intercept().  
            %  N.B.:  CI uses propagation of errors for quotients, CIs_for_quotient().

            if ~isempty(this.table_mito_brain_map_quotient_radar_)
                T = this.table_mito_brain_map_quotient_radar_;
                return
            end

            table_fqfn = fullfile(this.mito_profile_home, stackstr() + ".mat");
            if isfile(table_fqfn)
                ld = load(table_fqfn);
                T = ld.T;
                return
            end

            T_mbm = table(this);
            T_intercept = table_intercept(this);

            for mfp = this.mito_fps
                select = strcmp(T_mbm.ProfileLabel, mfp);
                ProfileLabel = T_mbm.ProfileLabel(select);
                T_mbm_1 = T_mbm(select, :);
                T_1 = this.CIs_for_quotient(T_mbm_1, T_intercept);
                T_1 = addvars(T_1, ProfileLabel, 'Before', 'Estimate');
                if strcmp(mfp, this.mito_fps(1))
                    T = T_1;
                else
                    T = [T; T_1];
                end
            end

            save(table_fqfn, "T");
            this.table_mito_brain_map_quotient_radar_ = T;
        end

        function s1 = call_mito_profiles(this, T, opts)
            arguments
                this mladni.MitoBrainMapRadar
                T table  %% with vars ProfileLabel, Estimate, CI
                opts.show logical = false
            end

            N_mito = length(this.mito_fps);
            e = zeros(N_mito, this.N_patterns);
            lci = zeros(N_mito, this.N_patterns);
            uci = zeros(N_mito, this.N_patterns);
            for im = 1:N_mito
                try
                    T_ = T(strcmp(T.ProfileLabel, this.mito_fps(im)), :);
                    if opts.show
                        disp(T_)
                    end

                    % Estimate
                    e(im,:) = T_.Estimate';
                    lci(im,:) = T_.LowerCI';
                    uci(im,:) = T_.UpperCI';

                catch ME
                    handwarning(ME)
                end
            end
        
            hfig = figure;
            amin = min(lci, [], "all");
            amax = max(uci, [], "all");
            fprintf("%s: P: amin->%g, amax->%g\n", stackstr(), amin, amax)
            c = this.suvr2viridis(min(e, [], 2));
            s1 = this.plot_groups_with_CI_bounds( ...
                e, lci, uci, ...
                AxesMin=amin-eps, AxesMax=amax+eps, ...
                AxesInterval=7, ...
                LineWidth=this.GROUP_LINEWIDTHS, ...
                MinorGrid="off", MinorGridInterval=[], ...
                legend=this.groupLabels, ...
                ti="", ...
                Color=c ...
            );
            saveFigure2(hfig, fullfile(this.figdir, stackstr()));
        end    

        function this = MitoBrainMapRadar(varargin)
            
            % Mito NIfTI
            this.mito_profile_home = fullfile( ...
                getenv("SINGULARITY_HOME"), ...
                "ADNI", "MitochondrialProfile", "2mm" ...
            );
            this.mito_fps = ["CI", "CII", "CIV", "MitoD", "MRC", "TRC"];
            this.mito_profiles = fullfile(this.mito_profile_home, this.mito_fps + ".nii.gz");

            % PoMC NIfTI
            this.basis_all = fullfile( ...
                getenv("SINGULARITY_HOME"), ...
                "ADNI", "NMF_FDG", "baseline_cn", "NumBases"+this.N_patterns, "OPNMF", "niiImg", "Basis_all.nii" ...
            );  % dim_4 == this.N_patterns == 24 by default

            % max of PoMC
            basis_ic = mlfourd.ImagingContext2(this.basis_all);
            this.basis_max = dipmax(basis_ic);

            % Basis_all.nii is sorted by FDG SUVR, so sorting indices are monotonic
            this.sorted_bases_ = 1:this.N_patterns;  % 1:24 by default

            % binary mask 
            this.intracranial_mask = mlfourd.ImagingContext2( ...
                fullfile(getenv("SINGULARITY_HOME"), "ADNI", "VolBin", "mask.nii.gz") ...
            );

            % misc
            this.figdir_ = this.mito_profile_home;
            this.groups = convertStringsToChars(this.mito_fps);
            this.groupLabels = convertStringsToChars(this.mito_fps);
            this.groups0 = convertStringsToChars(this.mito_fps);
            this.mergeDx = convertStringsToChars(this.mito_fps);
            this.matfile0 = fullfile(this.mito_profile_home, "MitoBrainMapRadar_table.mat");
            this.matfile_cohort = fullfile(this.mito_profile_home, "MitoBrainMapRadar_table.mat");
            this.matfile_cohort_sorted = fullfile(this.mito_profile_home, "MitoBrainMapRadar_table.mat");
            this.GROUP_LINEWIDTHS = 1.618*ones(1, length(this.mito_fps));
        end
    end

    %% PROTECTED

    properties (Access = protected)
        table_mito_brain_map_radar_
        table_mito_brain_map_quotient_radar_
    end

    methods (Access = protected)
        function s = plot_groups_with_CI_bounds(this, P, lowerCI, upperCI, opts)
            %% spider_plot_class_examples.m
            %  Example 5 with Excel-like radar charts and
            %  Example 9 with shaded areas
            %  lowerCI < P < upperCI

            arguments
                this mladni.MitoBrainMapRadar
                P double
                lowerCI double
                upperCI double
                opts.legend {mustBeText} = {''}
                opts.ti {mustBeTextScalar} = ''
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesPrecision double = 2
                opts.AxesScaling cell = {}
                opts.Color double = []
                opts.LineStyle {mustBeText} = '-'
                opts.LineWidth double = 1
                opts.AxesLabels = this.AxesLabelsNull
                opts.MinorGrid = "on"
                opts.MinorGridInterval double = [];
            end
            if isempty(opts.Color)
                opts.Color = this.suvr2viridis(P);
            end
            P = P(:, this.sorted_bases);
            lowerCI = lowerCI(:, this.sorted_bases);
            upperCI = upperCI(:, this.sorted_bases);
            axes_shaded_limits = { ...
                [lowerCI(1,:); upperCI(1,:)], ...
                [lowerCI(2,:); upperCI(2,:)], ...
                [lowerCI(3,:); upperCI(3,:)]};  % KLUDGE
            axes_shaded_colors = { ...
                opts.Color(1,:), opts.Color(2,:), opts.Color(3,:)};  % KLUDGE
            Nbases = size(P,2);
            
            % Delete variable in workspace if exists
            if exist('s', 'var')
                delete(s);
            end
            
            % Spider plot
            try
                s = spider_plot_class(P);
            catch ME
                handwarning(ME)
            end
            
            % Spider plot properties
            if ~isempty(opts.AxesMin) && ~isempty(opts.AxesMax)
                s.AxesLimits = [opts.AxesMin*ones(1,Nbases); opts.AxesMax*ones(1,Nbases)];
            end
            if ~isempty(opts.AxesInterval)
                s.AxesInterval = opts.AxesInterval;
            end
            s.MinorGrid = opts.MinorGrid;
            if ~isempty(opts.MinorGridInterval)
                s.MinorGridInterval = opts.MinorGridInterval;
            end
            s.AxesLabels = opts.AxesLabels;
            s.AxesShaded = 'on';
            s.AxesShadedLimits = axes_shaded_limits;
            s.AxesShadedColor = axes_shaded_colors;
            s.AxesShadedTransparency = this.ERR_ALPHA;
            s.AxesPrecision = opts.AxesPrecision;
            s.AxesDisplay = 'one';
            %s.FillOption = 'on';
            %s.FillTransparency = 0.1;
            s.Color = opts.Color;
            s.LineStyle = opts.LineStyle;
            s.LineWidth = opts.LineWidth;
            s.Marker = 'none';
            s.AxesFontSize = 11;
            s.LabelFontSize = 12;
            s.AxesColor = this.AXES_COLOR;
            s.AxesLabelsEdge = 'none';
            s.AxesLabelsRotate = 'on';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            if ~isemptytext(opts.legend)
                s.LegendLabels = opts.legend;
                s.LegendHandle.Location = 'northeastoutside';
                s.LegendHandle.FontSize = 13;
            end
            if ~isemptytext(opts.ti)
                title(opts.ti, 'FontSize', 14);
            end
        end
    end

    methods (Static, Access = protected)
        function [bootstrap_mean, ci_lower, ci_upper] = bootstrapROIMean(img_data, roi_mask, n_iterations, ci_level)
            % BOOTSTRAPROIMEAN Calculates bootstrapped mean and confidence intervals
            % for imaging data weighted by a probabilistic ROI mask
            %
            % Inputs:
            %   img_data      - 3D or 4D array containing imaging data
            %   roi_mask      - 3D array with values 0-1 representing ROI probabilities
            %   n_iterations  - Number of bootstrap iterations (default: 1000)
            %   ci_level      - Confidence interval level (default: 0.95)
            %
            % Outputs:
            %   bootstrap_mean - Mean of the bootstrapped distribution
            %   ci_lower       - Lower bound of confidence interval
            %   ci_upper       - Upper bound of confidence interval
            %
            % Claude 3.7 Sonnet for jjlee.wustl.edu@gmail.com.  2025 Apr 2.

            % Set defaults
            if nargin < 3
                n_iterations = 1000;
            end
            if nargin < 4
                ci_level = 0.95;
            end

            % Reshape data for easier handling
            img_size = size(img_data);
            if ndims(img_data) > 3
                % For 4D data (e.g., time series)
                n_voxels = prod(img_size(1:3));
                n_volumes = img_size(4);
                img_data = reshape(img_data, [n_voxels, n_volumes]);
                roi_mask = reshape(roi_mask, [n_voxels, 1]);
            else
                % For 3D data
                img_data = img_data(:);
                roi_mask = roi_mask(:);
                n_volumes = 1;
            end

            % Remove voxels with zero probability
            valid_idx = roi_mask > 0;
            img_data = img_data(valid_idx, :);
            roi_probs = roi_mask(valid_idx);

            % Normalize probabilities to sum to 1 (for sampling weights)
            roi_probs = roi_probs / sum(roi_probs);

            % Prepare for bootstrap
            n_valid_voxels = length(roi_probs);
            bootstrap_means = zeros(n_iterations, n_volumes);

            % Bootstrap sampling
            for iter = 1:n_iterations
                % Sample voxels with replacement according to ROI probabilities
                sample_idx = randsample(n_valid_voxels, n_valid_voxels, true, roi_probs);

                % Calculate weighted mean for this bootstrap sample
                sampled_data = img_data(sample_idx, :);
                sampled_weights = roi_probs(sample_idx);

                % Weighted mean across sampled voxels for each volume
                bootstrap_means(iter, :) = sum(sampled_data .* repmat(sampled_weights, 1, n_volumes), 1) / ...
                    sum(sampled_weights);
            end

            % Calculate statistics from bootstrap distribution
            alpha = 1 - ci_level;
            ci_lower_percentile = alpha/2 * 100;
            ci_upper_percentile = (1 - alpha/2) * 100;

            bootstrap_mean = mean(bootstrap_means, 1);
            ci_lower = prctile(bootstrap_means, ci_lower_percentile, 1);
            ci_upper = prctile(bootstrap_means, ci_upper_percentile, 1);

            % Reshape outputs to match original data dimensions if needed
            if n_volumes > 1
                bootstrap_mean = reshape(bootstrap_mean, [1, 1, 1, n_volumes]);
                ci_lower = reshape(ci_lower, [1, 1, 1, n_volumes]);
                ci_upper = reshape(ci_upper, [1, 1, 1, n_volumes]);
            end

        end
        function T = CIs_for_quotient(TA, TB)
            %  Let A be drawn from table TA, B be drawn from table TB.
            %  Let C = A/B.
            %  DeltaA ~ upper|lower confidence interval for A; equivalently for DeltaB and B.
            %  DeltaC/C = sqrt((DeltaA/A)^2 + (DeltaB/B)^2)
    
            assert(size(TA, 1) == size(TB, 1))
            N_rows = size(TA, 1);

            Estimate = [];
            LowerCI = [];
            UpperCI = [];
            for r = 1:N_rows
                e = TA.Estimate(r)/TB.Estimate(r);
                Estimate = [Estimate; e];
                
                lci = e * sqrt(TA.LowerCI(r)^2 / TA.Estimate(r)^2 + TB.LowerCI(r)^2 / TB.Estimate(r)^2);
                LowerCI = [LowerCI; lci];
                uci = e * sqrt(TA.UpperCI(r)^2 / TA.Estimate(r)^2 + TB.UpperCI(r)^2 / TB.Estimate(r)^2);
                UpperCI = [UpperCI; uci];
            end

            T = table(Estimate, LowerCI, UpperCI);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
