classdef Adni < mladni.DataCuration & handle
    %% ADNI is the top-level builder for ADNI-related objects and data.
    %  See also mladni.Oasis3
    %  
    %  Created 06-Aug-2023 15:38:27 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.

    properties (Dependent)
        data_home
        nmf_fdg_home
    end

    methods % GET
        function g = get.data_home(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI");
        end
        function g = get.nmf_fdg_home(this)
            g = fullfile(this.data_home, "NMF_FDG");
        end
    end

    methods
        function this = Adni()
            this.demogr = mladni.AdniDemographics();
            this.nmf = mladni.NMF( ...
                data_home=this.data_home, ...
                selectedNumBases=this.selectedNumBases);
            this.nmfc = mladni.NMFCovariates();
            this.nmfh = mladni.NMFHierarchies(data_home=this.data_home);
            this.nmfr = mladni.NMFRadar();
        end

        function anticlust_ad_repeat(this)
            %% prepares folders in NMF_FDG for reproducibility runs of NMF

            t = this.demogr.table_ad(true, 'table_fdg5');
            Filelist = strrep(t.Filelist, "/home/usr", "/scratch");
            t.Filelist = Filelist;
            csv_fqfn = fullfile(this.nmf_fdg_home, "anticlust_ad_repeat.csv");
            writetable(t, csv_fqfn, WriteVariableNames=true)

            R_fqfn = fullfile(this.nmf_fdg_home, "anticlust_ad_repeat.R");
            % mysystem(sprintf("Rscript %s", R_fqfn))

            mladni.Anticlust.prepare_folders_for_VolBin(data_home=this.data_home)
        end

        function anticlust_cn_repeat(this)
            %% prepares folders in NMF_FDG for reproducibility runs of NMF

            t = this.demogr.table_cn(true, 'table_fdg5');
            Filelist = strrep(t.Filelist, "/home/usr", "/scratch");
            t.Filelist = Filelist;
            csv_fqfn = fullfile(this.nmf_fdg_home, "anticlust_cn_repeat.csv");
            writetable(t, csv_fqfn, WriteVariableNames=true)

            R_fqfn = fullfile(this.nmf_fdg_home, "anticlust_cn_repeat.R");
            % mysystem(sprintf("Rscript %s", R_fqfn))

            mladni.Anticlust.prepare_folders_for_VolBin(data_home=this.data_home)
        end

        function call(this)
            % mladni.AdniBidsT1w ~ generates T1w in bids
            % mladni.AdniBidsFdg ~ generates FDG in bids
            % mladni.FDG ~ computes brain mask (dlicv), warps T1w & FDG to MNI, computes PVE, applies brain mask

            % generate table with filelist for NMF
            % call(this.demogr);

            % create summary images of other baseline cohorts
            % call2(this.demogr)
            
            % NMF on cluster
            % $SINGULARITY_HOME/ADNI/VolBin/submit_20230526.sh

            % create anticlust filelists with RStudio, then run on cluster
            % this.anticlust_cn_repeat();
            % $SINGULARITY_HOME/ADNI/VolBin/submit_anticlust_20230526.sh

            % write cache "X.mat", then calculate reconstruction errors
            % call(this.nmf)

            % calculate component-weighted averages
            % call2(this.nmf)

            % build argmax maps
            for b = 2:2:24
                nmfc = mladni.NMFCovariates(selectedNumBases=b);
                nmfc.table_covariates_1stscan();
            end
            this.nmfh.build_argmax_maps;

            % check completeness
            % fdg5 = this.demogr.table_fdg5;
            % assert(sum(cellfun(@isempty, fdg5.Filelist)) == 0, stackstr())

            % create images of means, var, std, median, iqr of baseline_cn, longitudinal_*
            subgroups = [ ...
                "baseline_cn", "longitudinal_cn", "longitudinal_preclinical", ...
                "longitudinal_cdr_0p5_apos", "longitudinal_cdr_gt_0p5_apos", ...
                "longitudinal_cdr_gt_0_aneg"];         
            for isg = 1:length(subgroups)
                try
                    % this.nmf.build_stats_imaging(inputDir=fullfile(this.nmf_fdg_home, subgroups(isg)));                    
                catch ME
                    handwarning(ME)
                end
            end
            % this.nmf.build_table_variances(subgroups=subgroups);

            % radar plots
            % plot_beta0_to_beta1(this.nmfr);
            % call_patt_weighted_fdg(this.nmfr);
            % call_groups(this.nmfr);
            % call_sex(this.nmfr);
            % call_apoe4(this.nmfr);

            % create surfaces
            % this.nmf.build_surfaces()

            % GAM(M4) regressions with RStudio

        end

        function h = paper_plot_rec_error(this, opts)
            arguments
                this mladni.Adni
                opts.fileprefix {mustBeTextScalar} = stackstr()
                opts.K double {mustBeInteger} = 20 % count of bases computed by VolBin/*
                opts.subgroup {mustBeTextScalar} = 'cn'
                opts.ylab {mustBeText} = 'ylab'
            end
        end

        function h = paper_plot_reproducibility(this, opts)
            %% plots evaluations of reproducibility, using metrics for selecting granular models.
            %  For each granular model:
            %  - raincloud representation of bootstraps of the similarity between combinations of patterns;
            %    similarity is the inner product of voxels of patterns from split halves of imaging samples, 
            %    adjusted for optimal matching by the Hungarian algorithm; split halves are outputs of 
            %    anti-clustering on all available imaging samples
            %  - median of bootstrapped similarities, also using median to marginalize comparator patterns
            %  - median of adjusted Rand Index from the hard clustering of spatial patterns generated by a model;
            %    hard clusters are the argmax of non-negative, real-valued spatially distributed functions 
            %    parameterized by the ordinal patterns of a model
            %  - median of fractional improvements to reconstruction errors, the residuals between linear combinations
            %    of patterns and the median of available imaging samples

            arguments
                this mladni.Adni
                opts.ARI {mustBeNumeric}
                opts.overlap {mustBeNumeric}
                opts.fileprefix {mustBeTextScalar} = stackstr()
                opts.K double {mustBeInteger} = 20 % count of bases computed by VolBin/*
                opts.subgroup {mustBeTextScalar} = 'cn'
                opts.ylab {mustBeText} = 'Split-sample reproducibility'
            end
            sortedBasisNum = 2:2:2*opts.K;
            sortedBasisNames = cellfun(@num2str, num2cell(2:2:2*opts.K), UniformOutput=false);

            rm_raincloud(sortedBasisNum, sortedBasisNames, ARI');
            xlabel('Number of patterns in model space', 'fontsize', 30)
            ylabel(opts.ylab, 'fontsize', 30)
            legend({'pattern overlaps', 'Adj. Rand Index'}, 'fontsize', 20)
            set(gca,'fontsize',20)
            set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])

            outputDir_ = fullfile(this.nmf_fdg_home, sprintf('baseline_%s', subgroup), 'results');
            saveas(gcf, fullfile(outputDir_, [opts.fileprefix, '.fig']));
            saveas(gcf, fullfile(outputDir_, [opts.fileprefix, '.png']));
            saveas(gcf, fullfile(outputDir_, [opts.fileprefix, '.svg']));
        end

        function t = paper_table_census(this)
            Groups = { ...
                'no. ADNI FDG scans'; ...
                'and T1w scans within 1y'; ...
                'and Pons-vermis reference within 1y'; ...
                'and visually quality controlled'; ...
                'and disjoint from OASIS3'; ...
                'and select 1st FDG'};
            n_rows = length(Groups);
            tbl = { ...
                'table_fdg1'; ...
                'table_fdg2'; ...
                'table_fdg4'; ...
                'table_fdg4'; ...
                'table_fdg5'; ...
                'table_fdg5'};
            assert(length(tbl) == length(Groups))
            viz_excl = zeros(n_rows, 6);
            viz_excl(4,:) = nan;
            viz_excl(4,1:2) = 1; % <- MANUALLY CONFIGURED
            assert(all(viz_excl(:,1) == sum(viz_excl(:,2:end), 2, 'omitnan')))

            All_available = zeros(n_rows, 1);
            CDR_0_aneg = zeros(n_rows, 1);
            CDR_0_apos = zeros(n_rows, 1);
            CDR_0p5_apos = zeros(n_rows, 1);
            CDR_gt_0p5_apos = zeros(n_rows, 1);
            CDR_gt_0_aneg = zeros(n_rows, 1);
            for r = 1:n_rows
                select_1st = r == n_rows;

                sz = size(this.demogr.table_all(select_1st, tbl{r}));
                All_available(r) = sz(1) - viz_excl(r,1); %#ok<*AGROW> 

                sz = size(this.demogr.table_cn(select_1st, tbl{r}));
                CDR_0_aneg(r) = sz(1) - viz_excl(r,2);

                sz = size(this.demogr.table_preclinical(select_1st, tbl{r}));
                CDR_0_apos(r) = sz(1) - viz_excl(r,3);

                sz = size(this.demogr.table_cdr_0p5_apos(select_1st, tbl{r}));
                CDR_0p5_apos(r) = sz(1) - viz_excl(r,4);

                sz = size(this.demogr.table_cdr_gt_0p5_apos(select_1st, tbl{r}));
                CDR_gt_0p5_apos(r) = sz(1) - viz_excl(r,5);

                sz = size(this.demogr.table_cdr_gt_0_aneg(select_1st, tbl{r}));
                CDR_gt_0_aneg(r) = sz(1) - viz_excl(r,6);
            end
            
            t = table(Groups, All_available, CDR_0_aneg, CDR_0_apos, CDR_0p5_apos, CDR_gt_0p5_apos, CDR_gt_0_aneg);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
