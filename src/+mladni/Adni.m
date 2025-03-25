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
                N_patterns=this.N_patterns, ...
                selected_spans=2:2:40);
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

            %% spatial normalizations

            % mladni.AdniBidsT1w ~ generates T1w in bids
            % mladni.AdniBidsFdg ~ generates FDG in bids
            % mladni.FDG ~ computes brain mask (dlicv), warps T1w & FDG to MNI, computes PVE, applies brain mask

            %% generate table with filelist for NMF
            
            % call(this.demogr);

            %% create summary images of other baseline cohorts
            
            % call2(this.demogr)
            
            %% NMF on cluster
            
            % $SINGULARITY_HOME/ADNI/VolBin/submit_20230526.sh

            %% create anticlust filelists with RStudio, then run on cluster
            
            % this.anticlust_cn_repeat();
            % $SINGULARITY_HOME/ADNI/VolBin/submit_anticlust_20230526.sh

            %% write cache "X.mat", calculate reconstruction errors, calculate component-weighted averages
            
            % call(this.nmf)

            % fdg5 = this.demogr.table_fdg5;
            % assert(sum(cellfun(@isempty, fdg5.Filelist)) == 0, stackstr())  % check completeness            

            %% build tables for diagnostic groups

            % for p = 2:2:40
            %     this.nmfc.N_patterns = p;
            %     this.nmfc.writetables();
            % end

            %% build argmax maps

            % this.nmfh.build_argmax_maps();
            % this.nmfh.build_table_for_ggalluvial2();

            %% run alluvial_across_ranks_tom2.Rmd in RStudio

            %% create images of means, var, std, median, iqr of baseline_cn, longitudinal_*

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
            % this.plot_model_selections()

            %% radar plots, unwrapped radar plots, preparations for chord plots

            % this.nmf.build_for_brainsmash();
            
            %% RStudio ADNI/NMF_FDG/baseline_cn/NumBases24/patterns_of_neurodegeneration_20240904.Rmd.qmd

            % call_intercept(this.nmfr);
            % call_groups(this.nmfr);
            % call_sex(this.nmfr);
            % call_apoe4(this.nmfr);
            % plot_beta0_to_beta1(this.nmfr);

            %% chord plots

            % ns = mladni.Neurosynth();
            % ns.build_stats_from_logs(tag="Neurosynth120");
            % ns.build_stats_from_logs(tag="Neurosynth27");
            
            %% run chords_for_neurosynth.qmd

            %% GPPM

            %% RStudio ADNI/NMF_FDG/baseline_cn/NumBases24/patterns_for_gppm.Rmd; 
            %% builds components/NMFCovariates_table_covariates_longitudinal_adjpatterns.csv

            % table_gppm_metarois(this.nmfc)
            % table_gppm_patterns(this.nmfc)
            % table_gppm_patterns_metaroi(this.nmfc)
            % table_gppm_patterns_1403(this.nmfc)
            
            %% run pycharm packages

            %% more QC

            %% create surfaces

            % this.nmf.build_surfaces()

            %% GAM(M4) regressions with RStudio




        end

        function call_hcp(this)
            %% sub-168S6085_ses-20171128111047_trc-FDG_proc-CASU-ponsvermis_orient-rpi_pet_on_T1w_Warped_dlicv.nii.gz  % youngest baseline
            %  
            %  https://ida.loni.usc.edu/pages/access/search.jsp


        end

        function h = plot_model_selections(this, opts)
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
                opts.subgroup {mustBeTextScalar} = 'cn'
            end

            removeColsWithNans = @mladni.NMF.removeColsWithNans;

            nmf_a = this.nmf;
            results_dir = fullfile(nmf_a.nmf_fdg_home, "baseline_" + opts.subgroup, "results");
            ld = load(fullfile(results_dir, "NMF_build_repeated_reproducibility.mat"));
            ARI_a = ld.ARI;
            overlap_numer_a = ld.overlap_numer;
            ld = load(fullfile(results_dir, "NMF_build_repeated_reproducibility_improvement_numer.mat"));
            improvement_numer_a = ld.improvement_numer;
            ld = load(fullfile(results_dir, "NMF_build_repeated_calc_rec_error.mat"));
            grad_numer_a = diff(cell2mat(ld.rec_errors), 1);

            nmf_o = mladni.NMF(data_home=fullfile(getenv("SINGULARITY_HOME"), "OASIS3"));
            results_dir = fullfile(nmf_o.nmf_fdg_home, "baseline_" + opts.subgroup, "results");
            ld = load(fullfile(results_dir, "NMF_build_repeated_reproducibility.mat"));
            ARI_o = ld.ARI;
            overlap_numer_o = ld.overlap_numer;
            ld = load(fullfile(results_dir, "NMF_build_repeated_reproducibility_improvement_numer.mat"));
            improvement_numer_o = ld.improvement_numer;
            ld = load(fullfile(results_dir, "NMF_build_repeated_calc_rec_error.mat"));
            grad_numer_o = diff(cell2mat(removeColsWithNans(ld.rec_errors)), 1);
            ld = load(fullfile(results_dir, "NMF_build_repeated_reproducibility3.mat"));
            ARI_ao = ld.ARI;
            overlap_numer_ao = ld.overlap_numer;

            plt = mladni.Plot();
            span_model = 2:2:40;
            d_span_model = 364.119 / 38;
            span_model_1 = (span_model - span_model(1)) * d_span_model;
            cred = cbrewer2('Reds', 9);
            cgreen = cbrewer2('YlGn', 9);
            cblue = cbrewer2('Blues', 9);
            crain = {cblue(6,:), cgreen(6,:), cred(6,:)};
            cpdf = {cblue(3,:), cgreen(3,:), cred(3,:)};
            cline = {cblue(8,:), cgreen(8,:), cred(8,:)};

            % %% ARIs
             
            figure(Position=nmf_a.fig_position);
            hold on;

            plt.rm_plot(span_model(1:size(ARI_a,1)), ARI_a', ...
                line_width=5, ...
                colours=cpdf{1}, ...
                line_colour=cline{1}, ...
                do_flip=false); 
            plt.rm_plot(span_model(1:size(ARI_o,1)), ARI_o', ...
                line_width=5, ...
                colours=cpdf{2}, ...
                line_colour=cline{2}, ...
                do_flip=false); 
            plt.rm_plot(span_model(1:size(ARI_ao,1)), ARI_ao', ...
                line_width=5, ...
                colours=cpdf{3}, ...
                line_colour=cline{3}, ...
                do_flip=false); 

            hold off;
            ylabel(["Reproducibility of hard clusters", "(adjusted Rand index)"]);
            xlabel("Number of patterns in model space");
            fontsize(scale=plt.fontsize_scale)
            legend(["ADNI IQR", "ADNI median", "OASIS3 IQR", "OASIS3 median", "ADNI-OASIS3 IQR", "ADNI-OASIS3 median"])

            %% inner products
             
            figure(Position=nmf_a.fig_position);
            hold on;

            data = {overlap_numer_a', overlap_numer_o', overlap_numer_ao'};
            plt.rm_raincloud(data, ...
                line_width=4, ...
                rain_colours=crain, ...
                pdf_colours=cpdf, ...
                line_colour=cline);

            hold off;
            xlabel(["Reproducibility of matched patterns", "(inner products)"]);
            ylabel("Number of patterns in model space");
            fontsize(scale=plt.fontsize_scale)
            legend(["ADNI PDF", "ADNI bootstraps", "OASIS3 PDF", "OASIS3 bootstraps", "ADNI-OASIS3 PDF", "ADNI-OASIS3 bootstraps"])

            %% reconstruction errors

            % figure(Position=nmf_a.fig_position);
            % hold on;
            % 
            % plt.rm_plot(span_model(1:size(improvement_numer_a,1)), improvement_numer_a', ...
            %     line_width=5, ...
            %     colours=c{1}(end - 1, :), ...
            %     line_colour=c{1}(end, :), ...
            %     do_flip=false); 
            % plt.rm_plot(span_model(1:size(improvement_numer_o,1)), improvement_numer_o', ...
            %     line_width=5, ...
            %     colours=c{2}(end - 1, :), ...
            %     line_colour=c{2}(end, :), ...
            %     do_flip=false); 
            % 
            % hold off;
            % ylabel("Cumulative improvement of reconstruction fidelity");
            % xlabel("Number of patterns in model space");
            % fontsize(scale=plt.fontsize_scale)
            % legend(["ADNI IQR", "ADNI median", "OASIS3 IQR", "OASIS3 median"])

            %% grad reconstruction errors

            % figure(Position=nmf_a.fig_position);
            % hold on;
            % 
            % plt.rm_plot(span_model(1:size(grad_numer_a,1)), grad_numer_a', ...
            %     line_width=5, ...
            %     colours=c{1}(end - 1, :), ...
            %     line_colour=c{1}(end, :), ...
            %     do_flip=false); 
            % plt.rm_plot(span_model(1:size(grad_numer_o,1)), grad_numer_o', ...
            %     line_width=5, ...
            %     colours=c{2}(end - 1, :), ...
            %     line_colour=c{2}(end, :), ...
            %     do_flip=false); 
            % 
            % hold off;
            % ylabel("Gradient of reconstruction error");
            % xlabel("Number of patterns in model space");
            % fontsize(scale=plt.fontsize_scale)
            % legend(["ADNI IQR", "ADNI median", "OASIS3 IQR", "OASIS3 median"])
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
