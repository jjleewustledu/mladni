classdef Adni < handle
    %% ADNI is the top-level builder for ADNI-related objects and data.
    %  See also mladni.Oasis3
    %  
    %  Created 06-Aug-2023 15:38:27 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        selectedNumBases = mladni.NMF.N_PATTERNS
    end

    properties
        ad
        nmf
        nmfc
        nmfr
    end

    properties (Dependent)
        workdir
    end

    methods % GET
        function g = get.workdir(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "NMF_FDG");
        end
    end

    methods
        function anticlust_cn_repeat(this)
            %% prepares folders in NMF_FDG for reproducibility runs of NMF

            t = this.ad.table_cn(true, 'table_fdg5');
            Filelist = strrep(t.Filelist, "/home/usr", "/scratch");
            t.Filelist = Filelist;
            csv_fqfn = fullfile(this.workdir, "anticlust_cn_repeat.csv");
            writetable(t, csv_fqfn, WriteVariableNames=true)

            R_fqfn = fullfile(this.workdir, "anticlust_cn_repeat.R");
            mysystem(sprintf("Rscript %s", R_fqfn))
        end
        function call(this)
            % mladni.AdniBidsT1w ~ generates T1w in bids
            % mladni.AdniBidsFdg ~ generates FDG in bids
            % mladni.FDG ~ computes brain mask (dlicv), warps T1w & FDG to MNI, computes PVE, applies brain mask

            % generate table with filelist for NMF
%            call(this.ad); 

            % create summary images of other baseline cohorts
%            call2(this.ad)
            
            % NMF on cluster
            % $SINGULARITY_HOME/ADNI/VolBin/submit_20230526.sh

            % create anticlust filelists with RStudio, then run on cluster
%            this.anticlust_cn_repeat();
            % $SINGULARITY_HOME/ADNI/VolBin/submit_anticlust_20230526.sh

            % write cache "X.mat", then calculate reconstruction errors
%            call(this.nmf)

            % estimate reproducibility
%            this.nmf.diagnose_reproducibility() % boxplots

            % calculate component-weighted averages
%            call2(this.nmf)
%            this.nmfc.table_covariates();
%            this.nmfc.table_covariates_1stscan();

            % check completeness
            fdg5 = this.ad.table_fdg5;
            assert(sum(@isempty, fdg5.Filelist) == 0, stackstr())

            % create images of means, var, std, median, iqr of baseline_cn, longitudinal_*
            subgroups = [ ...
                "baseline_cn", "longitudinal_cn", "longitudinal_preclinical", ...
                "longitudinal_cdr_0p5_apos", "longitudinal_cdr_gt_0p5_apos", ...
                "longitudinal_cdr_gt_0_aneg"];         
            for isg = 1:length(subgroups)
                try
                    this.nmf.build_stats_imaging(inputDir=fullfile(this.workdir, subgroups(isg)));                    
                catch ME
                    handwarning(ME)
                end
            end
            this.nmf.build_table_variances(subgroups);  

            % radar plots
            call_patt_weighted_fdg(this.nmfr);
            call_groups(this.nmfr);
            call_sex(this.nmfr);
            call_apoe4(this.nmfr);

            % create surfaces
%            this.nmf.build_surfaces()

            % GAM(M4) regressions with RStudio


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

                sz = size(this.ad.table_all(select_1st, tbl{r}));
                All_available(r) = sz(1) - viz_excl(r,1); %#ok<*AGROW> 

                sz = size(this.ad.table_cn(select_1st, tbl{r}));
                CDR_0_aneg(r) = sz(1) - viz_excl(r,2);

                sz = size(this.ad.table_preclinical(select_1st, tbl{r}));
                CDR_0_apos(r) = sz(1) - viz_excl(r,3);

                sz = size(this.ad.table_cdr_0p5_apos(select_1st, tbl{r}));
                CDR_0p5_apos(r) = sz(1) - viz_excl(r,4);

                sz = size(this.ad.table_cdr_gt_0p5_apos(select_1st, tbl{r}));
                CDR_gt_0p5_apos(r) = sz(1) - viz_excl(r,5);

                sz = size(this.ad.table_cdr_gt_0_aneg(select_1st, tbl{r}));
                CDR_gt_0_aneg(r) = sz(1) - viz_excl(r,6);
            end
            
            t = table(Groups, All_available, CDR_0_aneg, CDR_0_apos, CDR_0p5_apos, CDR_gt_0p5_apos, CDR_gt_0_aneg);
        end
        function this = Adni()
            this.ad = mladni.AdniDemographics();
            this.nmf = mladni.NMF(selectedNumBases=this.selectedNumBases);
            this.nmfc = mladni.NMFCovariates();
            this.nmfr = mladni.NMFRadar();
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
