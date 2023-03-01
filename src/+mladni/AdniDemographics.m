classdef AdniDemographics < handle
    %% queries "all_FDG_20211123_12_14_2021.csv"
    %
    %  14358 imaging items, all DCM, acquired 09/22/2005 - 10/13/2021
    %  1660 unique subjects, acquisition ages 55 - 96 years, 924 males, 736 females
    %  303 unique subjects progressed to AD
    %  visit numbers 2 - 101
    %  groups:      {'AD', 'CN', 'EMCI', 'LMCI', 'MCI', 'SMC'}, SMC ~ significant memory concern
    %  group scans: { 551   949   534     249     1341   111}
    %  
    %  3444 'Co-registered Dynamic'                              
    %  3444 'Co-registered, Averaged'                            
    %  3735 'Coreg, Avg, Standardized Image and Voxel Size'      
    %  3735 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution'}
    %
    % ImageDataID        Subject        Group      Sex     Age    Visit    Modality                          Description                              Type          AcqDate      Format     Downloaded
    % ____________    ______________    ______    _____    ___    _____    ________    _______________________________________________________    _____________    __________    _______    __________
    % 
    % {'I1480670'}    {'941_S_6962'}    {'AD'}    {'F'}    75       2      {'PET'}     {'Coreg, Avg, Standardized Image and Voxel Size'      }    {'Processed'}    08/10/2021    {'DCM'}    11/23/2021
    % {'I1480672'}    {'941_S_6962'}    {'AD'}    {'F'}    75       2      {'PET'}     {'Co-registered, Averaged'                            }    {'Processed'}    08/10/2021    {'DCM'}    11/23/2021
    % {'I1480671'}    {'941_S_6962'}    {'AD'}    {'F'}    75       2      {'PET'}     {'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution'}    {'Processed'}    08/10/2021    {'DCM'}    11/23/2021
    % {'I1480673'}    {'941_S_6962'}    {'AD'}    {'F'}    75       2      {'PET'}     {'Co-registered Dynamic'                              }    {'Processed'}    08/10/2021    {'DCM'}    11/23/2021
    % {'I1301874'}    {'941_S_6854'}    {'AD'}    {'M'}    86       2      {'PET'}     {'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution'}    {'Processed'}    02/24/2020    {'DCM'}    11/23/2021
    %  ...
    %
    %  Created 14-Dec-2021 13:51:20 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1809720 (R2021b) Update 1 for MACI64.  Copyright 2021 John J. Lee.
    
    methods (Static)
        function t = create_table_of_filenames(varargin)
            import mladni.AdniDemographics.glob_fdg_proc
            import mladni.AdniDemographics.glob_fdg_raw
            import mladni.AdniDemographics.glob_mpr

            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'subjectsDir', @isfolder)
            addParameter(ip, 'modality', '', @istext)
            addParameter(ip, 'table_filename', '', @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;
            if isempty(ipr.table_filename)
                ipr.table_filename = fullfile(ipr.subjectsDir, strcat(lower(ipr.modality), '_filenames.csv'));
            end
            
            switch lower(ipr.modality)
                case 'fdgproc'
                    c = glob_fdg_proc(ipr.subjectsDir);
                case 'fdgraw'
                    c = glob_fdg_raw(ipr.subjectsDir);
                case 't1'
                    c = glob_mpr(ipr.subjectsDir);
                otherwise
                    error('mladni:ValueError', 'AdniDemographics.create_table.ipr.modality->%s', ipr.modality)
            end
            t = cell2table(c, 'VariableNames', {'Filename'});
            writetable(t, ipr.table_filename);
        end
        function ic = fqfns2mean(fqfns, fileprefix)
            arguments
                fqfns cell
                fileprefix {mustBeTextScalar} = stackstr(2)
            end

            Nerr = 0;
            fqfns = fqfns(cellfun(@(x) isfile(x), fqfns));
            N = length(fqfns);
            ic = mlfourd.ImagingContext2(fqfns{1}) ./ N;
            for f = 2:N
                try
                ic = ic + mlfourd.ImagingContext2(fqfns{f}) ./ N;
                ic.fileprefix = fileprefix;
                catch 
                    Nerr = Nerr + 1;
                end
            end
            ic.fileprefix = strcat(ic.fileprefix, '_mean');
            ic.filepath = pwd;

            fprintf('%s: Nerr->%g\n', stackstr(), Nerr)
        end
        function g = glob_fdg_raw(varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'subjectsDir', @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            g = glob(fullfile(ipr.subjectsDir, '*_S_*', '*Raw*', '*', 'I*', ''));
            for ig = 1:length(g)
                hdrs = glob(fullfile(g{ig}, '*.i.hdr'));
                g{ig} = hdrs{1};
            end
            g = ascol(g);
        end
        function g = glob_fdg_proc(varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'subjectsDir', @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            g = glob(fullfile(ipr.subjectsDir, '*_S_*', 'Co*', '*', 'I*', ''));
            for ig = 1:length(g)
                dcms = glob(fullfile(g{ig}, '*.dcm'));
                g{ig} = dcms{1};
            end
            g = ascol(g);
        end
        function g = glob_mpr(varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'subjectsDir', @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            g = glob(fullfile(ipr.subjectsDir, '*_S_*', '*', '*', 'I*', ''));
            for ig = 1:length(g)
                niis = glob(fullfile(g{ig}, '*.nii*'));
                if ~isempty(niis)
                    g{ig} = niis{1};
                end
                dcms = glob(fullfile(g{ig}, '*.dcm'));
                if ~isempty(dcms)
                    g{ig} = dcms{1};
                end
            end
            g = ascol(g);
        end
        function g = table2fqfns(T, opts)
            % e.g., sub-*_ses-*_trc-FDG_proc-CASU*_orient-rpi_pet_on_T1w_Warped_dlicv.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_Warped.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_0_Warped.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_1_Warped.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_2_Warped.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_detJ.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_dlicv_detJ_Warped.nii.gz

            arguments
                T table = []
                opts.globbing {mustBeTextScalar} = 'sub-*_ses-*_trc-FDG_proc-CASU*_orient-rpi_pet_on_T1w_Warped_dlicv.nii.gz';
            end
            if isempty(T)
                ad = mladni.AdniDemographics;
                T = ad.table_fdg1();
            end

            g = {};
            for row = 1:size(T,1)
                try
                    sub_re = regexp(T{row, 'Subject'}{1}, '(?<pre>\d{3})_S_(?<rid>\d{4})', 'names');
                    sub = sprintf('sub-%sS%s', sub_re.pre, sub_re.rid);
                    ses_dt = T{row, 'AcqDate'};
                    ses_dt.Format = 'yyyyMMdd';
                    ses = strcat('ses-', char(ses_dt));
                    g_ = glob(fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', sub, ses, 'pet', opts.globbing));
                    if isempty(g_)
                        g_ = {''};
                    end
                    g = [g; g_]; %#ok<AGROW> 
                catch ME
                    g = [g; {''}];
                    handwarning(ME)
                end
            end
        end
    end

    properties
        description
        study_design
    end

    properties (Constant)
        AGLOBS = { ...
                  'sub-*_ses-*_trc-FDG_proc-CASU-ponsvermis_orient-rpi_pet_on_T1w_Warped_dlicv.nii.gz', ...
                  'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_dlicv_detJ_Warped.nii.gz', ...
                  'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_detJ.nii.gz', ...
                  'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_0_Warped.nii.gz', ...
                  'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_1_Warped.nii.gz', ...
                  'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_2_Warped.nii.gz', ...
                  'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_Warped.nii.gz'}
        categories = {'CN', 'EMCI', 'MCI', 'LMCI', 'SMC', 'AD'};
        LABELS = { ...
                  'sub-cn_ses-all_trc-FDG_pet_on_T1w_Warped_dlicv', ...
                  'sub-cn_ses-all_T1w_dlicv_detJ_Warped', ...
                  'sub-cn_ses-all_T1w_brain_detJ', ...
                  'sub-cn_ses-all_T1w_brain_pve_0_Warped', ...
                  'sub-cn_ses-all_T1w_brain_pve_1_Warped', ...
                  'sub-cn_ses-all_T1w_brain_pve_2_Warped', ...
                  'sub-cn_ses-all_T1w_brain_Warped'}
        LABELS_G = {'fdg' 'dlicv-detJ' 'brain-detJ' 'pve0' 'pve1' 'pve2' 'brain'}
        subgroups = { ...
                'cn', 'preclinical', ...
                'cdr_0p5_aneg_emci', 'cdr_0p5_aneg_lmci', ...
                'cdr_0p5_apos_emci', 'cdr_0p5_apos_lmci', ...
                'cdr_0p5_aneg', 'cdr_0p5_apos', 'cdr_0p5_anan', ...
                'cdr_gt_0p5_apos'}
        viscode2_months = [0 6 12 18 24 36 48 54 60 66 72 78 84 90 96 108 120 126 132 138 144 150];
    end

    properties (Dependent)
        adni_merge
        amyloid_file
        cdr_file
        fdg_orig_file
        fdg_proc1_file
        fdg_proc2_file
        fdg1_file
        fdgproc_filenames_file
        fqfns_scans_gt_1y %% N == 191
        pet_c3_file
        pet_meta_adni1_file
        pet_meta_adnigo2_file
        pet_meta_adni3_file
        pet_meta_list_file
        pet_qc_file
        subjects
        ucberkeleyfdg_file
    end

    methods % GET
        function g = get.adni_merge(this)
            g = this.adni_merge_;
        end
        function g = get.amyloid_file(~)
            g = fullfile(getenv("ADNI_HOME"), 'studydata', 'ucberkeley_av45_fbb_skinny.csv');
            % unique RID ~ 1325
        end
        function g = get.cdr_file(~)
            g = fullfile(getenv("ADNI_HOME"), "studydata", "CDR.csv"); %"cdr_20220602.csv");
            % unique RID ~ 3418
        end
        function g = get.fdg_orig_file(~)
            g = fullfile(getenv("ADNI_HOME"), "staging", "FDG_original_20211118", "all_FDG_20211118_11_18_2021.csv");
            % unique Subject ~ 1662
        end
        function g = get.fdg_proc1_file(~)
            g = fullfile(getenv("ADNI_HOME"), "staging", "FDG_processed_20211123", "all_FDG_20211123_12_14_2021.csv");
            % unique Subject ~ 1660, unique ImageDataID ~ 14358
        end
        function g = get.fdg_proc2_file(~)
            g = fullfile(getenv("ADNI_HOME"), "staging", "FDG_processed_20220104", "prev_failed_FDG_20220104_1_04_2022.csv");
            % unique Subject ~ 5, included in fdg_proc1_file
        end
        function g = get.fdg1_file(~)
            g = fullfile(getenv("ADNI_HOME"), "staging", "FDG_processed_20220104", "mladni_AdniDemographics_fdg1.csv");
        end
        function g = get.fdgproc_filenames_file(~)
            g = fullfile(getenv("ADNI_HOME"), "bids", "rawdata", "rosters", "fdgproc_filenames.csv");
        end      
        function g = get.fqfns_scans_gt_1y(this)
            if isempty(this.fqfns_scans_gt_1y_)
                ld = load(fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', 'fqfns_scans_gt_1y.mat'));
                this.fqfns_scans_gt_1y_ = ld.fqfns_scans_gt_1y; 
            end
            g = this.fqfns_scans_gt_1y_;
        end
        function g = get.pet_c3_file(~)
            g = fullfile(getenv("ADNI_HOME"), "studydata", "PETC3.csv");
            % unique RID ~ 395, unique LONIUID ~ 407
        end
        function g = get.pet_meta_adni1_file(~)
            g = fullfile(getenv("ADNI_HOME"), "studydata", "PETMETA_ADNI1.csv");
            % unique RID ~ 420
        end
        function g = get.pet_meta_adnigo2_file(~)
            g = fullfile(getenv("ADNI_HOME"), "studydata", "PETMETA_ADNIGO2.csv");
            % unique RID ~ 1212
        end
        function g = get.pet_meta_adni3_file(~)
            g = fullfile(getenv("ADNI_HOME"), "studydata", "PETMETA3.csv");
            % unique RID ~ 641
        end
        function g = get.pet_meta_list_file(~)
            g = fullfile(getenv("ADNI_HOME"), "studydata", "PET_META_LIST.csv");
            % unique Subject ~ 2028, unique ImageID ~ 50843
        end
        function g = get.pet_qc_file(~)
            g = fullfile(getenv("ADNI_HOME"), "studydata", "PETQC.csv");
            % unique RID ~ 1413, unique LONIUID ~ 3950
        end
        function g = get.subjects(this)
            if isempty(this.subjects_)
                tbl_fdg = this.table_fdg();
                this.subjects_ = unique(tbl_fdg.Subject);
            end
            g = this.subjects_;
        end
        function g = get.ucberkeleyfdg_file(~)
            g = fullfile(getenv("ADNI_HOME"), "studydata", "UCBERKELEYFDG_03_25_22.csv");
        end
    end

    methods
        function this = AdniDemographics(desc, design)
            %% ADNIDEMOGRAPHICS
            %  Args:
            %      desc {mustBeTextScalar} = 'Co-registered Dynamic' | 'Coreg, Avg, Standardized Image and Voxel Size'
            %      design {mustBeTextScalar} = 'cross-sectional' | 'longitudinal'

            arguments
                desc {mustBeTextScalar} = 'Coreg, Avg, Standardized Image and Voxel Size'
                design {mustBeTextScalar} = 'cross-sectional'
            end

            this.description = desc;
            this.study_design = design;
            warning('off', 'MATLAB:table:ModifiedAndSavedVarnames') 
            this.adni_merge_ = mladni.AdniMerge();
        end

        function t = reduce_cross_tables(this, t1, t2, opts)
            %% Reduce cross-sectional tables

            arguments
                this mladni.AdniDemographics  
                t1 table
                t2 table
                opts.vars cell = {'RID', 'EXAMDATE', 'TRACER', 'AmyloidStatus', 'SUMMARYSUVR_WHOLECEREBNORM', 'SUMMARYSUVR_COMPOSITE_REFNORM'} % {pivot_var, ..., kept_var, ...}
                opts.pivot_vars {mustBeText} = {'RID', 'EXAMDATE'}
            end

            t_ = reduce_long_tables(this, t1, t2, vars=opts.vars, pivot_vars=opts.pivot_vars);
            t  = table();
            urid = unique(t_.RID);
            for idx = 1:length(urid)
                trid_ = t_(t_.RID == urid(idx),:);
                trid_ = sortrows(trid_, 'EXAMDATE');
                t = [t; trid_(1,:)]; %#ok<AGROW> 
            end
        end
        function t = reduce_long_tables(this, t1, t2, opts)
            %% Reduce longitudinal tables

            arguments
                this mladni.AdniDemographics %#ok<INUSA> 
                t1 table
                t2 table
                opts.vars cell = {'RID', 'EXAMDATE', 'TRACER', 'AmyloidStatus', 'SUMMARYSUVR_WHOLECEREBNORM', 'SUMMARYSUVR_COMPOSITE_REFNORM'} % {pivot_var, ..., kept_var, ...}
                opts.pivot_vars {mustBeText} = {'RID', 'EXAMDATE'}
            end

            not_vars_t1 = t1.Properties.VariableNames(~contains(t1.Properties.VariableNames, opts.vars));
            for nv = not_vars_t1
                t1.(nv{1}) = [];
            end
            not_vars_t2 = t2.Properties.VariableNames(~contains(t2.Properties.VariableNames, opts.vars));
            for nv = not_vars_t2
                t2.(nv{1}) = [];
            end
            t = [t1; t2];
            t = sortrows(t, opts.pivot_vars);
            
        end
        function t = table_amyloid(this, varargin)
            %% https://adni.bitbucket.io/reference/docs/UCBERKELEYAV45/UCBERKELEY_AV45_Methods_04.25.2022.pdf
            %  https://adni.bitbucket.io/reference/docs/UCBERKELEYFBB/UCBerkeley_FBB_Methods_04.25.2022.pdf

            if ~isempty(this.amyloid_)
                t = this.amyloid_;
            end
            if isempty(this.amyloid_) && isfile(this.amyloid_file)
                t = readtable(this.amyloid_file);
                t.EXAMDATE = datetime(t.EXAMDATE);
                this.amyloid_ = t;
            end
            if isempty(this.amyloid_) && ~isfile(this.amyloid_file)
                % read tables of av45, florbetaben
                pth = fileparts(this.amyloid_file);                
                av45 = fullfile(pth, 'UCBERKELEYAV45_04_26_22.csv');
                assert(isfile(av45))
                t_av45 = readtable(av45);
                t_av45.TRACER = repmat("av45", [size(t_av45,1) 1]);
                fbb = fullfile(pth, 'UCBERKELEYFBB_04_26_22.csv');
                assert(isfile(fbb))
                t_fbb = readtable(fbb);
                t_fbb.TRACER = repmat("fbb", [size(t_fbb,1) 1]);

                % reduce tables according to study design
                switch this.study_design
                    case 'cross-sectional'
                        t_av45.AmyloidStatus = double(t_av45{:, 'SUMMARYSUVR_WHOLECEREBNORM'} >= 1.11);
                        t_fbb.AmyloidStatus = double(t_fbb{:, 'SUMMARYSUVR_WHOLECEREBNORM'} >= 1.08);
                        t = this.reduce_cross_tables(t_av45, t_fbb);
                    case 'longitudinal'
                        t_av45.AmyloidStatus = double(t_av45{:, 'SUMMARYSUVR_COMPOSITE_REFNORM'} >= 0.78);
                        t_fbb.AmyloidStatus = double(t_fbb{:, 'SUMMARYSUVR_COMPOSITE_REFNORM'} >= 0.74);
                        t = this.reduce_long_tables(t_av45, t_fbb);
                    otherwise
                        error('mladni:ValueError', 'this.study_design->%s', this.study_design)
                end
                writetable(t, this.amyloid_file);
            end
            
            % slice the table with varargin
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_cdr(this, varargin)
            if isempty(this.cdr_)
                this.cdr_ = this.cdrRevisions(readtable(this.cdr_file));
            end
            t = this.cdr_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_damyloid(this, varargin)
            t = table;
            t_ = this.table_fdg1;
            t_ = t_(contains(t_.Description, 'Uniform')',:); % CASU

            select = cellfun(@(x) ~isempty(x), t_.MergeDxBl); % select availabe MergeDxBl
            t_ = t_(select,:);

            t_.damyloid = nan(size(t_.AmyloidStatus));

            % t <- Subject's last MergeExamDate, damyloid
            isub = 0;
            for S = unique(t_.Subject)'
                select = strcmp(S{1}, t_.Subject);
                t_S = t_(select, :);
                t_S_ascend = sortrows(t_S, 'MergeExamDate', 'ascend');
                if t_S_ascend{1, 'AmyloidStatus'} == 0 
                    t_S_ascend{1, 'damyloid'} = 0;
                    if size(t_S_ascend, 1) > 1
                        t_S_ascend{2:end, 'damyloid'} = t_S_ascend{2:end, 'AmyloidStatus'} > 0;
                    end
                end
                if t_S_ascend{1, 'AmyloidStatus'} == 1 
                    t_S_ascend{1, 'damyloid'} = 0;
                    if size(t_S_ascend, 1) > 1
                        t_S_ascend{2:end, 'damyloid'} = 0;
                    end
                end
                [~,imax] = max(t_S_ascend.MergeExamDate);
                isub = isub + 1;
                t(isub, :) = t_S_ascend(imax, :);
            end
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_dCDGLOBAL(this, varargin)
            t = table;
            t_ = this.table_fdg1;
            t_ = t_(contains(t_.Description, 'Uniform')',:); % CASU

            select = cellfun(@(x) ~isempty(x), t_.MergeDxBl); % select availabe MergeDxBl
            t_ = t_(select,:);

            t_.dCDGLOBAL = nan(size(t_.CDGLOBAL));

            % t <- Subject's last MergeExamDate, dCDGLOBAL
            isub = 0;
            for S = unique(t_.Subject)'
                select = strcmp(S{1}, t_.Subject);
                t_S = t_(select, :);
                t_S_ascend = sortrows(t_S, 'MergeExamDate', 'ascend');

                t_S_ascend{1, 'dCDGLOBAL'} = 0;
                if size(t_S_ascend, 1) > 1
                    t_S_ascend{2:end, 'dCDGLOBAL'} = t_S_ascend{2:end, 'CDGLOBAL'} - t_S_ascend{1, 'CDGLOBAL'};
                end

                [~,imax] = max(t_S_ascend.dCDGLOBAL);
                isub = isub + 1;
                t(isub, :) = t_S_ascend(imax, :);
            end
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_dMergeCdrsb(this, varargin)
            t = table;
            t_ = this.table_fdg1;
            t_ = t_(contains(t_.Description, 'Uniform')',:); % CASU

            select = cellfun(@(x) ~isempty(x), t_.MergeDxBl); % select availabe MergeDxBl
            t_ = t_(select,:);

            t_.dMergeCdrsb = nan(size(t_.MergeCdrsb));

            % t <- Subject's last MergeExamDate, dMergeCdrsb
            isub = 0;
            for S = unique(t_.Subject)'
                select = strcmp(S{1}, t_.Subject);
                t_S = t_(select, :);
                t_S_ascend = sortrows(t_S, 'MergeExamDate', 'ascend');

                t_S_ascend{1, 'dMergeCdrsb'} = 0;
                if size(t_S_ascend, 1) > 1
                    t_S_ascend{2:end, 'dMergeCdrsb'} = t_S_ascend{2:end, 'MergeCdrsb'} - t_S_ascend{1, 'MergeCdrsb'};
                end

                [~,imax] = max(t_S_ascend.dMergeCdrsb);
                isub = isub + 1;
                t(isub, :) = t_S_ascend(imax, :);
            end
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_duration(this, varargin)
            t = table;
            t_ = this.table_fdg1;
            t_ = t_(contains(t_.Description, 'Uniform')',:); % CASU

            select = cellfun(@(x) ~isempty(x), t_.MergeDxBl); % select availabe MergeDxBl
            t_ = t_(select,:);

            Duration = years(t_.MergeExamDate - t_.MergeExamDateBl);
            t_.Duration = Duration;

            % t <- Subject's max Duration
            isub = 0;
            for S = unique(t_.Subject)'
                select = strcmp(S{1}, t_.Subject);
                t_S = t_(select, :);
                [~,imax] = max(t_S.Duration);
                isub = isub + 1;
                t(isub, :) = t_S(imax, :);
            end
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_fdg(this, varargin)
            if isempty(this.fdg_)
                this.fdg_ = readtable(this.fdg_proc1_file);
            end
            t = this.fdg_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_fdg_select_description(this, varargin)
            t = this.fdg_(strcmp(this.fdg_.Description, this.description), :);
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end  
        function t = table_fdg1(this, varargin)

            % cached in memory
            if ~isempty(this.fdg1_)
                t = this.fdg1_;
                if ~isempty(varargin)
                    t = t(varargin{:});
                end       
                return
            end

            % cached on filesystem
            cache_file = fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', strcat(clientname(true, 2), '.mat'));
            if isfile(cache_file)
                ld = load(cache_file);
                this.fdg1_ = ld.table_fdg1;
                t = this.fdg1_;
                if ~isempty(varargin)
                    t = t(varargin{:});
                end      
                return
            end
            
            % build caches
            this.fdg1_ = this.buildTableFdg1();
            table_fdg1 = this.fdg1_;
            save(cache_file, 'table_fdg1');
            t = this.fdg1_;            
            if ~isempty(varargin)
                t = t(varargin{:});
            end            
        end
        function t = table_fdgproc_filenames(this, varargin)
            if isempty(this.fdgproc_filenames_)
                this.fdgproc_filenames_ = readtable(this.fdgproc_filenames_file);
            end
            t = this.fdgproc_filenames_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_firstscan(this, varargin)
            %% lazy init with table_fdg()

            if isempty(this.firstscan_)
                t_ = this.table_fdg1();
                t = t_(strcmp(t_.Subject, this.subjects{1}), :); 
                t = t(t.AcqDate == min(t.AcqDate), :); 
                for ti = 2:length(this.subjects)
                    u = t_(strcmp(t_.Subject, this.subjects{ti}), :); % pick subject
                    v = u(u.AcqDate == min(u.AcqDate), :); % pick first scan                
                    t = [t; v]; %#ok<AGROW> % append 
                end
                this.firstscan_ = t;
            else
                t = this.firstscan_;
            end
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_lastscan(this, varargin)
            %% lazy init with table_fdg()

            if isempty(this.lastscan_)
                t_ = this.table_fdg1();
                t = t_(strcmp(t_.Subject, this.subjects_{1}), :); 
                t = t(t.AcqDate == max(t.AcqDate), :); 
                for ti = 2:length(this.subjects_)
                    u = t_(strcmp(t_.Subject, this.subjects_{ti}), :); % pick subject
                    v = u(u.AcqDate == max(u.AcqDate), :); % pick first scan                
                    t = [t; v]; %#ok<AGROW> % append 
                end
                this.lastscan_ = t;
            else
                t = this.lastscan_;
            end
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end        
        function t = table_pet_c3(this, varargin)
            if isempty(this.pet_c3_)
                this.pet_c3_ = readtable(this.pet_c3_file);
            end
            t = this.pet_c3_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_pet_qc(this, varargin)
            if isempty(this.pet_qc_)
                this.pet_qc_ = readtable(this.pet_qc_file);
            end
            t = this.pet_qc_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_select_available(~, t_in, varname)
            avail = cellfun(@(x) ~isempty(x), t_in.(varname));
            t = t_in(avail, :);
        end
        function t = table_select_description(this, t_in, varargin)
            t = t_in(strcmp(t_in.Description, this.description), :);
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end 
        function t = table_ucberkeleyfdg(this, varargin)
            if isempty(this.ucberkeleyfdg_)
                t_ = readtable(this.ucberkeleyfdg_file);            
                t_.EXAMDATE = datetime(t_.EXAMDATE);
                this.ucberkeleyfdg_ = t_;
            end
            t = this.ucberkeleyfdg_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end        
        function     write_imageids(~, tbl, fname)
            %% Creates a file with comma-separated imageIDs. 
            %  For LONI IDA workflows, see also:
            %  https://github.com/sotiraslab/adni_collections/blob/main/collections/allt1/directions.ipynb

            c = tbl.ImageDataID;
            if isnumeric(tbl.ImageDataID{1})
                c = num2cell(c);
                c = cellfun(@(x) num2str(x), c, 'UniformOutput', false);
            end
            s = strjoin(c, ',');
            fid = fopen(fname, 'w');
            fprintf(fid, s); 
            fclose(fid);
        end

        %% subgroups defined
        %  1660 baseline FDG scans with CASU

        function     call_subgroups_ic_means(this)
            %% sequentially calls create_subgroups(), create_ic_means()

            f1 = @this.create_subgroups;
            aglobs = this.AGLOBS;
            lblg = this.LABELS_G;
            f2 = @this.create_ic_means;
            lbl = this.LABELS;
            for idx = 1:length(aglobs)
%                f1(aglobs{idx}, lblg{idx});
                f2(lbl{idx}, lblg{idx});
            end
        end
        function     create_ic_means(this, label, labelg)
            %% 
            %  Args:
            %     label {mustBeTextScalar} = 'sub-cn_ses-all_T1w_dlicv_detJ_Warped'
            %     label {mustBeTextScalar} = 'sub-cn_ses-all_T1w_brain_detJ'
            %     label {mustBeTextScalar} = 'sub-cn_ses-all_T1w_brain_pve_0_Warped'
            %     label {mustBeTextScalar} = 'sub-cn_ses-all_T1w_brain_Warped'
            %     label {mustBeTextScalar} = 'sub-cn_ses-all_trc-FDG_pet_on_T1w_Warped_dlicv'

            arguments
                this mladni.AdniDemographics
                label {mustBeTextScalar} = 'sub-cn_ses-all_trc-FDG_pet_on_T1w_Warped_dlicv'
                labelg {mustBeTextScalar} = 'fdg'
            end
            
            load(strcat(labelg, '_fqfns.mat')) %#ok<*LOAD> 
            
            for sg = this.subgroups

                % save diagnostic subgroups
                ic.(sg{1}) = this.fqfns2mean(fqfns.(sg{1}));
                ic.(sg{1}).fileprefix = strrep(label, this.subgroups{1}, sg{1});
                ic.(sg{1}).save();   
            end

            for sg = this.subgroups(2:end)

                % differences, view_qc()
                icd.(sg{1}) = ic.(sg{1}) - ic.(this.subgroups{1});
                icd.(sg{1}).fileprefix = strrep(label, this.subgroups{1}, strcat('D', sg{1}));
                %ic.(this.subgroups{1}).view_qc(icd.(sg{1}));
                ic.(this.subgroups{1}).save_qc(icd.(sg{1}));
            end
        end 
        function     create_subgroups(this, aglob, label)
            %% 
            %  Args:
            %     aglob {mustBeTextScalar} = 'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_dlicv_detJ_Warped.nii.gz'
            %     aglob {mustBeTextScalar} = 'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_detJ.nii.gz'
            %     aglob {mustBeTextScalar} = 'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_0_Warped.nii.gz'
            %     aglob {mustBeTextScalar} = 'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_Warped.nii.gz'
            %     aglob {mustBeTextScalar} = 'sub-*_ses-*_trc-FDG_proc-CASU-ponsvermis_orient-rpi_pet_on_T1w_Warped_dlicv.nii.gz'
            %     label {mustBeTextScalar} = 'dlicv_detJ'
            %     label {mustBeTextScalar} = 'brain_detJ'
            %     label {mustBeTextScalar} = 'pve0'
            %     label {mustBeTextScalar} = 'brain'
            %     label {mustBeTextScalar} = 'fdg'

            arguments
                this mladni.AdniDemographics
                aglob {mustBeTextScalar} = 'sub-*_ses-*_trc-FDG_proc-CASU-ponsvermis_orient-rpi_pet_on_T1w_Warped_dlicv.nii.gz'
                label {mustBeTextScalar} = 'fdg'
            end
            
            to_rem = mybasename(this.fqfns_scans_gt_1y);
            
            for sg = this.subgroups
                table_ = strcat('table_', sg{1});
                T.(sg{1}) = this.(table_); % call this.table_*()
                fqfns.(sg{1}) = this.table2fqfns(T.(sg{1}), globbing=aglob); % call this.table2fqfns()
                
                Filename = fqfns.(sg{1}); % rename for addvars()
                T.(sg{1}) = addvars(T.(sg{1}), Filename, 'Before', 1);
                T.(sg{1}) = renamevars(T.(sg{1}), "Var1", "Filename");
                
                fqfn_ = fqfns.(sg{1});
                fqfns.(sg{1})  = fqfn_(~contains(fqfn_, to_rem) & ~isempty(fqfn_)); % remove to_rem, scans > 1y separated
                T_ = T.(sg{1});
                T.(sg{1}) = T_(~contains(T_.Filename, to_rem) & ~isempty(T_.Filename), :);

                writetable(T.(sg{1}), sprintf('table_%s.csv', sg{1}), WriteVariableNames=true)
            end

            save(sprintf('%s_fqfns.mat', label), 'fqfns')
            save('T.mat', 'T')
        end
        function t = table_cn(this)
            %% N = 248

            first = this.table_firstscan();
            first = this.table_select_description(first);
            t = first(first.AmyloidStatus == 0 & first.CDGLOBAL == 0, :);
        end
        function t = table_preclinical(this)
            %% N = 133

            first = this.table_firstscan();
            first = this.table_select_description(first);
            t = first(first.AmyloidStatus == 1 & first.CDGLOBAL == 0, :);
        end
        function t = table_cdr_0p5_aneg(this)
            %% N = 272

            first = this.table_firstscan();
            first = this.table_select_description(first);
            t = first(first.CDGLOBAL == 0.5 & first.AmyloidStatus == 0, :);
        end
        function t = table_cdr_0p5_aneg_emci(this)
            %% N = 170

            first = this.table_firstscan();
            first = this.table_select_description(first);
            first = this.table_select_available(first, 'MergeDxBl');
            t = first(first.CDGLOBAL == 0.5 & first.AmyloidStatus == 0 & strcmp(first.MergeDxBl, 'EMCI'), :);
        end
        function t = table_cdr_0p5_aneg_lmci(this)
            %% N = 86

            first = this.table_firstscan();
            first = this.table_select_description(first);
            first = this.table_select_available(first, 'MergeDxBl');
            t = first(first.CDGLOBAL == 0.5 & first.AmyloidStatus == 0 & strcmp(first.MergeDxBl, 'LMCI'), :);
        end
        function t = table_cdr_0p5_apos(this)
            %% N = 385

            first = this.table_firstscan();
            first = this.table_select_description(first);
            t = first(first.CDGLOBAL == 0.5 & first.AmyloidStatus == 1, :);
        end
        function t = table_cdr_0p5_apos_emci(this)
            %% N = 152

            first = this.table_firstscan();
            first = this.table_select_description(first);
            first = this.table_select_available(first, 'MergeDxBl');
            t = first(first.CDGLOBAL == 0.5 & first.AmyloidStatus == 1 & strcmp(first.MergeDxBl, 'EMCI'), :);
        end
        function t = table_cdr_0p5_apos_lmci(this)
            %% N = 170

            first = this.table_firstscan();
            first = this.table_select_description(first);
            first = this.table_select_available(first, 'MergeDxBl');
            t = first(first.CDGLOBAL == 0.5 & first.AmyloidStatus == 1 & strcmp(first.MergeDxBl, 'LMCI'), :);
        end
        function t = table_cdr_0p5_anan(this)
            %% N == 190

            first = this.table_firstscan();
            first = this.table_select_description(first);
            t = first(first.CDGLOBAL == 0.5 & isnan(first.AmyloidStatus), :);
        end
        function t = table_cdr_gt_0p5_apos(this)
            %% N == 91

            first = this.table_firstscan();
            first = this.table_select_description(first);
            t = first(first.CDGLOBAL > 0.5 & first.AmyloidStatus == 1, :);
        end
        function t = table_cdr_gt_0p5_aneg(this)
            %% N = 14

            first = this.table_firstscan();
            first = this.table_select_description(first);
            t = first(first.CDGLOBAL > 0.5 & first.AmyloidStatus == 0, :);
        end

        %% visualizations

        function h = histogram_acqdate(this)
            figure
            h = histogram(this.table_fdg_select_description.AcqDate, "Normalization","count");
            xlabel("acquisition date")
            ylabel("count")
        end
        function h = histogram_age(this)
            figure
            h = histogram(this.table_fdg_select_description.Age, "Normalization","count");
            xlabel("age at acquisition")
            ylabel("count")
        end
        function h = histogram_duration_enrollment(this)
            figure
            last = this.table_lastscan(this.table_fdg1);
            first = this.table_firstscan(this.table_fdg1);
            assert(all(strcmp(last.Subject, first.Subject)))
            dt = calmonths(between(first.AcqDate, last.AcqDate));
            h = histogram(dt, "Normalization","count","BinWidth",6);
            xlabel("duration of enrollment (months)")
            ylabel("unique subjects")

            figure
            plot(sort(dt))
            xlabel("unique subject")
            ylabel("duration of enrollment (months)")
        end
        function h = histogram_group(this)
            figure
            t = this.table_fdg_select_description();
            c = categorical(t.Group, this.categories, 'Ordinal', true);
            summary(c)
            histogram(c)
            xlabel("clinical group")
            ylabel("count")
        end
        function h = histogram_progression(this, group)
            %% histogram including only subjects who progressed through a specified group.
            %  Args:
            %      group (text): drawn from this.categories

            assert(istext(group))

            t_ = this.table_fdg_select_description();
            t_group = t_(strcmp(t_.Group, group), :);
            subjects_group = unique(t_group.Subject);

            t = t_(strcmp(t_.Subject, subjects_group{1}), :);
            for ti = 2:length(subjects_group)
                u = t_(strcmp(t_.Subject, subjects_group{ti}), :); % pick subject
                t = [t; u]; %#ok<AGROW> % append
            end

            figure
            c = categorical(t.Group, this.categories, 'Ordinal', true);
            summary(c)
            histogram(c)
            xlabel("clinical group for subjects who progressed to " + group)
            ylabel("count")
        end
        function h = histogram_visits(this)
            figure
            h = histogram(this.table_fdg_select_description.Visit, "Normalization","count","NumBins",101);
            xlabel("visit number")
            ylabel("count")
        end
        function h = boxplotVISCODE2(this)
            dur = calmonths(between(this.fdg1_.MergeExamDateBl, this.fdg1_.MergeExamDate));
            vc = categorical(this.fdg1_.VISCODE2, ...
                {'bl' 'm06' 'm12' 'm18' 'm24' 'm36' 'm48' 'm54' 'm60' 'm66' 'm72' 'm78' 'm84' 'm90' 'm96' 'm108' 'm120' 'm126' 'm132' 'm138' 'm144' 'm150'}, ...
                'Ordinal', true);
            figure
            boxplot(dur, vc)
            xlabel('VISCODE2')
            ylabel('duration of enrollment (months)')
            set(gca, 'FontSize', 16)
        end
        function h = csplotVISCODE2(this)
            dur = calmonths(between(this.fdg1_.MergeExamDateBl, this.fdg1_.MergeExamDate));
            vc = this.fdg1_.VISCODE2;
            vc(strcmp(vc, 'm06')) = {'m006'};
            vc(strcmp(vc, 'm12')) = {'m012'};
            vc(strcmp(vc, 'm18')) = {'m018'};
            vc(strcmp(vc, 'm24')) = {'m024'};
            vc(strcmp(vc, 'm36')) = {'m036'};
            vc(strcmp(vc, 'm48')) = {'m048'};
            vc(strcmp(vc, 'm54')) = {'m054'};
            vc(strcmp(vc, 'm60')) = {'m060'};
            vc(strcmp(vc, 'm66')) = {'m066'};
            vc(strcmp(vc, 'm72')) = {'m072'};
            vc(strcmp(vc, 'm78')) = {'m078'};
            vc(strcmp(vc, 'm84')) = {'m084'};
            vc(strcmp(vc, 'm90')) = {'m090'};
            vc(strcmp(vc, 'm96')) = {'m096'};
            vc(strcmp(vc, 'm108')) = {'m108'};
            vc(strcmp(vc, 'm120')) = {'m120'};
            vc(strcmp(vc, 'm126')) = {'m126'};
            vc(strcmp(vc, 'm132')) = {'m132'};
            vc(strcmp(vc, 'm138')) = {'m138'};
            vc(strcmp(vc, 'm144')) = {'m144'};
            vc(strcmp(vc, 'm150')) = {'m150'};
            %labels = {'bl' 'm06' 'm12' 'm18' 'm24' 'm36' 'm48' 'm54' 'm60' 'm66' 'm72' 'm78' 'm84' 'm90' ...
            %          'm96' 'm108' 'm120' 'm126' 'm132' 'm138' 'm144' 'm150'};
            figure
            CategoricalScatterplot(dur, vc, ...
                'BoxColor', [0.8471 0.8627 0.8392])
            xlabel('VISCODE2')
            ylabel('duration of enrollment (months)')
            set(gca, 'FontSize', 16)
        end
        function h = csplotVISCODE2_for_fdg(this)
            dur = calmonths(between(this.fdg1_.MergeExamDateBl, this.fdg1_.AcqDate));
            vc = this.fdg1_.VISCODE2;
            vc(strcmp(vc, 'm06')) = {'m006'};
            vc(strcmp(vc, 'm12')) = {'m012'};
            vc(strcmp(vc, 'm18')) = {'m018'};
            vc(strcmp(vc, 'm24')) = {'m024'};
            vc(strcmp(vc, 'm36')) = {'m036'};
            vc(strcmp(vc, 'm48')) = {'m048'};
            vc(strcmp(vc, 'm54')) = {'m054'};
            vc(strcmp(vc, 'm60')) = {'m060'};
            vc(strcmp(vc, 'm66')) = {'m066'};
            vc(strcmp(vc, 'm72')) = {'m072'};
            vc(strcmp(vc, 'm78')) = {'m078'};
            vc(strcmp(vc, 'm84')) = {'m084'};
            vc(strcmp(vc, 'm90')) = {'m090'};
            vc(strcmp(vc, 'm96')) = {'m096'};
            vc(strcmp(vc, 'm108')) = {'m108'};
            vc(strcmp(vc, 'm120')) = {'m120'};
            vc(strcmp(vc, 'm126')) = {'m126'};
            vc(strcmp(vc, 'm132')) = {'m132'};
            vc(strcmp(vc, 'm138')) = {'m138'};
            vc(strcmp(vc, 'm144')) = {'m144'};
            vc(strcmp(vc, 'm150')) = {'m150'};
            %labels = {'bl' 'm06' 'm12' 'm18' 'm24' 'm36' 'm48' 'm54' 'm60' 'm66' 'm72' 'm78' 'm84' 'm90' ...
            %          'm96' 'm108' 'm120' 'm126' 'm132' 'm138' 'm144' 'm150'};
            figure
            CategoricalScatterplot(dur, vc, ...
                'BoxColor', [0.8471 0.8627 0.8392])
            xlabel('VISCODE2')
            ylabel('time of FDG from baseline (months)')
            set(gca, 'FontSize', 16)
        end
        function c = viscode2(this)
            if isempty(this.fdg1_)
                this.table_fdg1();
            end
            c = unique(this.fdg1_.VISCODE2);
            % 'bl m06 m12 m18 m24 m36 m48 m54 m60 m66 m72 m78 m84 m90 m96 m108 m120 m126 m132 m138 m144 m150'
        end
    end

    %% PROTECTED

    properties (Access = protected)
        adni_merge_
        amyloid_
        cdr_
        fdg_
        fdg1_
        fdgproc_filenames_
        firstscan_
        fqfns_scans_gt_1y_
        lastscan_
        pet_c3_
        pet_qc_
        subjects_
        ucberkeleyfdg_
    end

    methods (Access = protected)
        function t = buildTableFdg1(this)
            
            t = this.table_fdg();
            sz = size(t.AcqDate);

            t.AcqDate = datetime(t.AcqDate);
            t.AmyloidStatus = nan(sz);

            t.MergeVisCode = cell(sz);
            t.MergeExamDate = NaT(sz);
            t.MergeDxBl = cell(sz);
            t.MergeAge = nan(sz);
            t.MergePtGender = cell(sz);
            t.MergePtEducat = nan(sz);
            t.MergePtEthCat = cell(sz);
            t.MergePtRacCat = cell(sz);
            t.MergePtMarry = cell(sz);
            t.MergeApoE4 = nan(sz);
            t.MergeFdg = nan(sz);
            t.MergePib = nan(sz);
            t.MergeAv45 = nan(sz);
            t.MergeAbeta = nan(sz);
            t.MergeTau = nan(sz);
            t.MergePTau = nan(sz);
            t.MergeCdrsb = nan(sz);
            t.MergeMmse = nan(sz);
            t.MergeDx = cell(sz);

            t.MergeExamDateBl = NaT(sz);
            t.MergeCdrsbBl = nan(sz);
            t.MergeAdas11Bl = nan(sz);
            t.MergeAdas13Bl = nan(sz);
            t.MergeAdasq4Bl = nan(sz);
            t.MergeMmseBl = nan(sz);
            t.MergeRavltImmediateBl = nan(sz);
            t.MergeRavltLearningBl = nan(sz);
            t.MergeRavltForgettingBl = nan(sz);
            t.MergeRavltPercForgettingBl = nan(sz);
            t.MergeEcogPtMemBl = nan(sz);
            t.MergeEcogPtLangBl = nan(sz);
            t.MergeEcogPtVisspatBl = nan(sz);
            t.MergeEcogPtPlanBl = nan(sz);
            t.MergeEcogPtOrganBl = nan(sz);
            t.MergeEcogPtDivattBl = nan(sz);
            t.MergeEcogPtTotalBl = nan(sz);
            t.MergeEcogSPMemBl = nan(sz);
            t.MergeEcogSPLangBl = nan(sz);
            t.MergeEcogSPVisspatBl = nan(sz);
            t.MergeEcogSPPlanBl = nan(sz);
            t.MergeEcogSPOrganBl = nan(sz);
            t.MergeEcogSPDivattBl = nan(sz);
            t.MergeEcogSPTotalBl = nan(sz);
            t.MergeAbetaBl = nan(sz);
            t.MergeTauBl = nan(sz);
            t.MergePTauBl = nan(sz);
            t.MergeFdgBl = nan(sz);
            t.MergePibBl = nan(sz);
            t.MergeAv45Bl = nan(sz);
            t.MergeYearsBl = nan(sz);
            t.MergeMonthBl = nan(sz);

            t.PonsVermis = nan(sz);
            t.VISCODE2 = cell(sz);
            t.USERDATE = NaT(sz);
            t.EXAMDATE = NaT(sz);
            t.CDMEMORY = nan(sz);
            t.CDORIENT = nan(sz);
            t.CDJUDGE = nan(sz);
            t.CDCOMMUN = nan(sz);
            t.CDHOME = nan(sz);
            t.CDCARE = nan(sz);
            t.CDGLOBAL = nan(sz);
            t.Phase = cell(sz);
            t.ID = nan(sz);
            t.RID = nan(sz);
            t.SITEID = nan(sz);

            t = movevars(t, 'RID', 'After', 'Subject');

            for s = asrow(this.subjects)
                try
                    t_s1 = t(strcmp(t.Subject, s{1}), :);
                    merge_s1 = this.adni_merge.table_merge(strcmp(this.adni_merge.table_merge.PTID, s{1}), ':');
                    if isempty(merge_s1)
                        continue
                    end
                    rid = merge_s1.RID(1);
                    assert(all(rid == merge_s1.RID))
                    cdr_s1 = this.table_cdr(this.table_cdr.RID == rid, ':');
                    try
                        amyloid_s1 = this.table_amyloid(this.table_amyloid.RID == rid, ':');
                        amy_sorted = sort(amyloid_s1.EXAMDATE);
                        amy_initial = amy_sorted(1) - years(1); % init av45 - 1 year
                        amy_future = t_s1.AcqDate > amy_initial; % future fdg following av45
                    catch
                    end
                    berkeley_s1 = this.table_ucberkeleyfdg(this.table_ucberkeleyfdg.RID == rid, ':');
                    berkeley_s1 = berkeley_s1(contains(berkeley_s1.ROINAME, 'pons-vermis'), :);
    
                    numrows_s1 = size(t_s1, 1);
                    for acqdi = 1:numrows_s1
                        acqdate = t_s1.AcqDate(acqdi);

                        [~,merge_near] = min(abs(acqdate - merge_s1.EXAMDATE));
                        t_s1(acqdi, 'MergeVisCode') = merge_s1(merge_near, 'VISCODE');
                        t_s1(acqdi, 'MergeExamDate') = merge_s1(merge_near, 'EXAMDATE'); % ------- EXAMDATE -------
                        t_s1(acqdi, 'MergeDxBl') = merge_s1(merge_near, 'DX_bl');
                        t_s1(acqdi, 'MergeAge') = merge_s1(merge_near, 'AGE');
                        t_s1(acqdi, 'MergePtGender') = merge_s1(merge_near, 'PTGENDER');
                        t_s1(acqdi, 'MergePtEducat') = merge_s1(merge_near, 'PTEDUCAT');
                        t_s1(acqdi, 'MergePtEthCat') = merge_s1(merge_near, 'PTETHCAT');
                        t_s1(acqdi, 'MergePtRacCat') = merge_s1(merge_near, 'PTRACCAT');
                        t_s1(acqdi, 'MergePtMarry') = merge_s1(merge_near, 'PTMARRY');
                        t_s1(acqdi, 'MergeApoE4') = merge_s1(merge_near, 'APOE4');
                        t_s1(acqdi, 'MergeFdg') = merge_s1(merge_near, 'FDG');
                        t_s1(acqdi, 'MergePib') = merge_s1(merge_near, 'PIB');
                        t_s1(acqdi, 'MergeAv45') = merge_s1(merge_near, 'AV45');
                        t_s1(acqdi, 'MergeAbeta') = merge_s1(merge_near, 'ABETA');
                        t_s1(acqdi, 'MergeTau') = merge_s1(merge_near, 'TAU');
                        t_s1(acqdi, 'MergePTau') = merge_s1(merge_near, 'PTAU');
                        t_s1(acqdi, 'MergeCdrsb') = merge_s1(merge_near, 'CDRSB');
                        t_s1(acqdi, 'MergeMmse') = merge_s1(merge_near, 'MMSE');
                        t_s1(acqdi, 'MergeDx') = merge_s1(merge_near, 'DX');

                        t_s1(acqdi, 'MergeExamDateBl') = merge_s1(merge_near, 'EXAMDATE_bl');
                        t_s1(acqdi, 'MergeCdrsbBl') = merge_s1(merge_near, 'CDRSB_bl');
                        t_s1(acqdi, 'MergeAdas11Bl') = merge_s1(merge_near, 'ADAS11_bl');
                        t_s1(acqdi, 'MergeAdas13Bl') = merge_s1(merge_near, 'ADAS13_bl');
                        t_s1(acqdi, 'MergeAdasq4Bl') = merge_s1(merge_near, 'ADASQ4_bl');
                        t_s1(acqdi, 'MergeMmseBl') = merge_s1(merge_near, 'MMSE_bl');
                        t_s1(acqdi, 'MergeRavltImmediateBl') = merge_s1(merge_near, 'RAVLT_immediate_bl');
                        t_s1(acqdi, 'MergeRavltLearningBl') = merge_s1(merge_near, 'RAVLT_learning_bl');
                        t_s1(acqdi, 'MergeRavltForgettingBl') = merge_s1(merge_near, 'RAVLT_forgetting_bl');
                        t_s1(acqdi, 'MergeRavltPercForgettingBl') = merge_s1(merge_near, 'RAVLT_perc_forgetting_bl');
                        t_s1(acqdi, 'MergeEcogPtMemBl') = merge_s1(merge_near, 'EcogPtMem_bl');
                        t_s1(acqdi, 'MergeEcogPtLangBl') = merge_s1(merge_near, 'EcogPtLang_bl');
                        t_s1(acqdi, 'MergeEcogPtVisspatBl') = merge_s1(merge_near, 'EcogPtVisspat_bl');
                        t_s1(acqdi, 'MergeEcogPtPlanBl') = merge_s1(merge_near, 'EcogPtPlan_bl');
                        t_s1(acqdi, 'MergeEcogPtOrganBl') = merge_s1(merge_near, 'EcogPtOrgan_bl');
                        t_s1(acqdi, 'MergeEcogPtDivattBl') = merge_s1(merge_near, 'EcogPtDivatt_bl');
                        t_s1(acqdi, 'MergeEcogPtTotalBl') = merge_s1(merge_near, 'EcogPtTotal_bl');
                        t_s1(acqdi, 'MergeEcogSPMemBl') = merge_s1(merge_near, 'EcogSPMem_bl');
                        t_s1(acqdi, 'MergeEcogSPLangBl') = merge_s1(merge_near, 'EcogSPLang_bl');
                        t_s1(acqdi, 'MergeEcogSPVisspatBl') = merge_s1(merge_near, 'EcogSPVisspat_bl');
                        t_s1(acqdi, 'MergeEcogSPPlanBl') = merge_s1(merge_near, 'EcogSPPlan_bl');
                        t_s1(acqdi, 'MergeEcogSPOrganBl') = merge_s1(merge_near, 'EcogSPOrgan_bl');
                        t_s1(acqdi, 'MergeEcogSPDivattBl') = merge_s1(merge_near, 'EcogSPDivatt_bl');
                        t_s1(acqdi, 'MergeEcogSPTotalBl') = merge_s1(merge_near, 'EcogSPTotal_bl');
                        t_s1(acqdi, 'MergeAbetaBl') = merge_s1(merge_near, 'ABETA_bl');
                        t_s1(acqdi, 'MergeTauBl') = merge_s1(merge_near, 'TAU_bl');
                        t_s1(acqdi, 'MergePTauBl') = merge_s1(merge_near, 'PTAU_bl');
                        t_s1(acqdi, 'MergeFdgBl') = merge_s1(merge_near, 'FDG_bl');
                        t_s1(acqdi, 'MergePibBl') = merge_s1(merge_near, 'PIB_bl');
                        t_s1(acqdi, 'MergeAv45Bl') = merge_s1(merge_near, 'AV45_bl');
                        t_s1(acqdi, 'MergeYearsBl') = merge_s1(merge_near, 'Years_bl');
                        t_s1(acqdi, 'MergeMonthBl') = merge_s1(merge_near, 'Month_bl');

                        [~,cdr_near] = min(abs(acqdate - cdr_s1.USERDATE));
                        t_s1(acqdi, 'VISCODE2') = cdr_s1{cdr_near, 'VISCODE2'};
                        t_s1{acqdi, 'USERDATE'} = cdr_s1{cdr_near, 'USERDATE'};
                        t_s1{acqdi, 'EXAMDATE'} = cdr_s1{cdr_near, 'EXAMDATE'}; % ------- EXAMDATE -------
                        t_s1(acqdi, 'CDMEMORY') = cdr_s1(cdr_near, 'CDMEMORY');
                        t_s1(acqdi, 'CDORIENT') = cdr_s1(cdr_near, 'CDORIENT');
                        t_s1(acqdi, 'CDJUDGE') = cdr_s1(cdr_near, 'CDJUDGE');
                        t_s1(acqdi, 'CDCOMMUN') = cdr_s1(cdr_near, 'CDCOMMUN');
                        t_s1(acqdi, 'CDHOME') = cdr_s1(cdr_near, 'CDHOME');
                        t_s1(acqdi, 'CDCARE') = cdr_s1(cdr_near, 'CDCARE');
                        t_s1(acqdi, 'CDGLOBAL') = cdr_s1(cdr_near, 'CDGLOBAL');
                        t_s1(acqdi, 'Phase') = cdr_s1(cdr_near, 'Phase');
                        t_s1(acqdi, 'ID') = cdr_s1(cdr_near, 'ID');
                        t_s1(acqdi, 'RID') = cdr_s1(cdr_near, 'RID');
                        t_s1(acqdi, 'SITEID') = cdr_s1(cdr_near, 'SITEID');   
                        
                        if isnat(t_s1{acqdi, 'EXAMDATE'}) && ~isnat(acqdate) % ------- EXAMDATE -------
                            t_s1{acqdi, 'EXAMDATE'} = acqdate;
                        end
                        
                        % amyloid status 
                        try
                            [~,amyloid_near] = min(abs(acqdate - amyloid_s1.EXAMDATE));
                            t_s1(acqdi, 'AmyloidStatus') = amyloid_s1(amyloid_near, 'AmyloidStatus');
                            t_s1(acqdi, 'SUMMARYSUVR_WHOLECEREBNORM') = amyloid_s1(amyloid_near, 'SUMMARYSUVR_WHOLECEREBNORM');
                            t_s1(acqdi, 'SUMMARYSUVR_COMPOSITE_REFNORM') = amyloid_s1(amyloid_near, 'SUMMARYSUVR_COMPOSITE_REFNORM');
                        catch
                        end

                        % find pons-vermis reference values
                        try
                            [~,berkeley_near] = min(abs(acqdate - berkeley_s1.EXAMDATE));
                            t_s1(acqdi, 'PonsVermis') = berkeley_s1(berkeley_near, 'MEAN');
                        catch 
                        end
                    end
                    
                    t(strcmp(t.Subject, s{1}), :) = t_s1;
                catch ME
                    handwarning(ME)
                end
            end
            
            % update t with table_registry.EXAMDATE (cf. ADNI Google Group)
            t_reg = this.adni_merge.table_registry();
            for ti = 1:size(t, 1)
                t_ti = t(ti, :);
                ED1 = t_reg(t_reg.ID == t_ti.ID, :).EXAMDATE;
                if isscalar(ED1) && ~isnat(ED1)
                    t(ti, :).EXAMDATE = ED1;
                end
            end
        end
        function t = cdrRevisions(~, t)
            return
            
            v2 = t.VISCODE2;
            v2(strcmp(v2, 'sc')) = {'bl'};
            v2(strcmp(v2, 'f')) = {'bl'};
            t.VISCODE2 = v2;
            t.USERDATE = datetime(t.USERDATE);
            t.USERDATE2 = datetime(t.USERDATE2);
            t.EXAMDATE = datetime(t.EXAMDATE);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
