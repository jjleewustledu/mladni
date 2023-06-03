classdef AdniDemographics < handle
    %% Supports ADNI demographic data and other meta-data at granularity larger than single scan objects.
    %  
    %  Queries "all_FDG_20211123_12_14_2021.csv", "ADNIMERGE.csv", mladni.AdniMerge
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
    %  See also:  https://groups.google.com/g/adni-data
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
            ic = ic ./ (N - Nerr);
            ic.fileprefix = strcat(ic.fileprefix, '_mean');
            ic.filepath = pwd;

            fprintf('%s: Nerr->%g\n', stackstr(), Nerr)
        end
        function list_acqdate = fqfns2AcqDate(list_fqfns)
            assert(iscell(list_fqfns), stackstr)
            list_acqdate = cell(size(list_fqfns));
            for idx = 1:length(list_fqfns)
                re = regexp(list_fqfns{idx}, '\S+sub-(?<sub>\d{3}S\d{4})_ses-(?<acqdate>\d{14})\S+', 'names');
                list_acqdate{idx} = datetime(re.acqdate, InputFormat='yyyyMMddHHmmss');
            end
        end
        function list_sub = fqfns2Subject(list_fqfns)
            assert(iscell(list_fqfns), stackstr)
            list_sub = cell(size(list_fqfns));
            for idx = 1:length(list_fqfns)
                re = regexp(list_fqfns{idx}, '\S+sub-(?<sub>\d{3}S\d{4})_ses-\S+', 'names');
                list_sub{idx} = strrep(re.sub, 'S', '_S_');
            end
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
            % e.g., opts.globbing ->
            %       sub-*_ses-*_trc-FDG_proc-CASU*_orient-rpi_pet_on_T1w_Warped_dlicv.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_Warped.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_0_Warped.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_1_Warped.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_pve_2_Warped.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_detJ.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_dlicv_detJ_Warped.nii.gz

            %       sub-*_ses-*_trc-FDG_proc-CASU_orient-rpi_pet.nii.gz
            %       sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w.nii.gz

            arguments
                T table = []
                opts.globbing {mustBeTextScalar} = 'sub-*_ses-*_trc-FDG_proc-CASU_orient-rpi_pet_on_T1w.nii.gz';
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
        function t = table_paren(varargin)
            t = mladni.AdniMerge.table_paren(varargin{:});
        end
        function t = table_excluding_table(t1, t2)
            %% from t1, exclude all rows with Subject & AcqDate found anywhere in t2

            assert(any(contains(t1.Properties.VariableNames, 'Subject')), stackstr)
            assert(any(contains(t1.Properties.VariableNames, 'AcqDate')), stackstr)
            assert(any(contains(t2.Properties.VariableNames, 'Subject')), stackstr)
            assert(any(contains(t2.Properties.VariableNames, 'AcqDate')), stackstr)

            select = false(size(t1,1), 1);
            for row2 = 1:size(t2, 1)
                selectSubject = contains(t1.Subject, t2.Subject{row2});
                selectAcqDate = t1.AcqDate == t2.AcqDate(row2);
                select = select | (selectSubject & selectAcqDate);
            end
            t = t1(~select, :);
        end
    end

    properties
        datetime_separation_tol
        description
        filt_bl_first % select earliest AcqDate before any other selection rules
        study_design
        reuse_cache
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
        LABELS_G = {'fdg'} % {'dlicv-detJ' 'brain-detJ' 'pve0' 'pve1' 'pve2' 'brain'}
        subgroups = { ...
                'cn', ...
                'preclinical', ...
                'cdr_0p5_apos', ...
                'cdr_gt_0p5_apos', ...
                'cdr_ge_0p5_aneg'}
        viscode2_months = [0 6 12 18 24 36 48 54 60 66 72 78 84 90 96 108 120 126 132 138 144 150];
    end

    properties (Dependent)
        adni_merge
        description_acro
        fdg1_file
        fdg_orig_file
        fdg_proc1_file
        fdg_proc2_file
        fqfns_scans_gt_1y %% N == 191, FDG scans with T1w > 1y of separation
        merge
        subjects
        workdir
    end

    methods % GET
        function g = get.adni_merge(this)
            g = this.merge;
        end  
        function g = get.description_acro(this)
            switch this.description
                case 'Co-registered Dynamic'
                    g = 'CD';
                case 'Co-registered, Averaged'
                    g = 'CA';
                case 'Coreg, Avg, Standardized Image and Voxel Size'
                    g = 'CAS';
                case 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution'
                    g = 'CASU';
                otherwise
                    error('mladni:ValueError', stackstr())
            end
        end
        function g = get.fdg1_file(this)
            g = fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', ...
                sprintf('AdniDemographics_table_fdg1_%s.mat', this.study_design));
        end
        function g = get.fdg_orig_file(~)
            %% unique Subject ~ 1662

            g = fullfile(getenv("ADNI_HOME"), "staging", "FDG_original_20211118", "all_FDG_20211118_11_18_2021.csv");            
        end
        function g = get.fdg_proc1_file(~)
            %% unique Subject ~ 1660, unique ImageDataID ~ 14358

            g = fullfile(getenv("ADNI_HOME"), "staging", "FDG_processed_20211123", "all_FDG_20211123_12_14_2021.csv");            
        end
        function g = get.fdg_proc2_file(~)
            %% unique Subject ~ 5, included in fdg_proc1_file

            g = fullfile(getenv("ADNI_HOME"), "staging", "FDG_processed_20220104", "prev_failed_FDG_20220104_1_04_2022.csv");            
        end 
        function g = get.fqfns_scans_gt_1y(this)
            if isempty(this.fqfns_scans_gt_1y_)
                this.fqfns_scans_gt_1y_ = {'disabling_fqfns_scans_gt_1y'};
                % ld = load(fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', 'fqfns_scans_gt_1y.mat'));
                % this.fqfns_scans_gt_1y_ = ld.fqfns_scans_gt_1y; 
            end
            g = this.fqfns_scans_gt_1y_;
        end
        function g = get.merge(this)
            g = this.merge_;
        end
        function g = get.subjects(this)
            if isempty(this.subjects_)
                tbl_fdg = this.table_fdg();
                this.subjects_ = unique(tbl_fdg.Subject);
            end
            g = this.subjects_;
        end
        function g = get.workdir(~)
            g = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
        end
    end

    methods
        function this = AdniDemographics(opts)
            %% ADNIDEMOGRAPHICS
            %  Args:
            %      opts.desc {mustBeTextScalar} = 'Co-registered Dynamic' | 'Coreg, Avg, Standardized Image and Voxel Size'
            %      opts.design {mustBeTextScalar} = 'cross-sectional' | 'longitudinal'
            %      opts.reuse_cache logical = false
            %      opts.filt_bl_first logical = true
            %      opts.datetime_separation_tol duration = years(100)
            %      opts.propagate_amyloid_positivity logical = true

            arguments
                opts.desc {mustBeTextScalar} = 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution'
                opts.design {mustBeTextScalar} = 'longitudinal'
                opts.reuse_cache logical = true
                opts.filt_bl_first logical = false
                opts.datetime_separation_tol duration = years(1)
            end
            this.description = opts.desc;
            this.study_design = opts.design;
            this.reuse_cache = opts.reuse_cache;
            this.filt_bl_first = opts.filt_bl_first;   
            this.datetime_separation_tol = opts.datetime_separation_tol;
            warning('off', 'MATLAB:table:ModifiedAndSavedVarnames') 
            this.merge_ = mladni.AdniMerge(pwd, this.study_design, reuse_cache=opts.reuse_cache);
        end

        function T = addvars_separation(this, opts)
            %% separation := datetime(measurement) - datetime(fdg acqdate), in days
            %  must read all NIfTI in ADNI/bids/derivatives

            arguments
                this mladni.AdniDemographics
                opts.T table = this.table_fdg1
                opts.measurement {mustBeTextScalar} = 'T1w'
            end
            measurement = convertStringsToChars(opts.measurement);
            switch lower(measurement)
                case 't1w'
                    Nrows = size(opts.T, 1);
                    t1w_acqdate = NaT(Nrows, 1);
                    for row = 1:Nrows
                        try
                            T_row = opts.T(row,:);
                            fqfns = this.table2fqfns(T_row, ...
                                globbing='sub-*_ses-*_trc-FDG_proc-CASU_orient-rpi_pet_on_T1w.json');
                            j = jsondecode(fileread(fqfns{1}));
                            re = regexp(j.T1wFilename, '\S+/sub-(?<sub>\d{3}S\d{4})_ses-(?<acqdate>\d{14})\S+', 'names');
                            t1w_acqdate(row) = datetime(re.acqdate, InputFormat='yyyyMMddHHmmss');
                        catch ME
                            fprintf('skipping %s %s; %s\n', T_row.Subject{1}, T_row.AcqDate, ME.message);
                        end
                    end
                    separation = days(t1w_acqdate - opts.T.AcqDate);
                case 'cdglobal'
                    separation = days(opts.T.CDEXAMDATE - opts.T.AcqDate);
                case 'amyloid'
                    separation = days(opts.T.AmyloidExamDate - opts.T.AcqDate);
                otherwise
                    error('mladni:ValueError', '%s.measurement->%s', stackstr, measurement)
            end
            T = addvars(opts.T, separation, ...
                NewVariableNames={sprintf('SeparationOf%s%s', upper(measurement(1)), lower(measurement(2:end)))});
        end    
        function as = AmyloidStatus(this, acqdate, amyloid_s1)
            arguments
                this mladni.AdniDemographics
                acqdate datetime
                amyloid_s1 table % derived from this.table_amyloid, for single subject
            end

            assert(1 == length(unique(amyloid_s1.RID)), stackstr)
            amyloid_s1 = amyloid_s1(~isnat(amyloid_s1.EXAMDATE), :);
            amyloid_s1 = sortrows(amyloid_s1, 'EXAMDATE'); % ascending

            % no known amyloidosis through end of amy scanning
            if all(amyloid_s1.AmyloidStatus == 0)
                if acqdate <= amyloid_s1.EXAMDATE(end) + duration(this.datetime_separation_tol)
                    as = 0;
                    return
                else
                    as = nan;
                    return
                end
            end

            % known amyloidosis since start of amy scanning
            if amyloid_s1.AmyloidStatus(1)
                if acqdate >= amyloid_s1.EXAMDATE(1) - duration(this.datetime_separation_tol)
                    as = 1;
                    return
                else
                    as = nan;
                    return
                end
            end

            % at least one amy+ scan after first amy scan
            [~,idx] = max(amyloid_s1.AmyloidStatus); 
            as = double(acqdate >= amyloid_s1.EXAMDATE(idx) - duration(this.datetime_separation_tol));
        end
        function as = AmyloidStatusCS(this, acqdate, amyloid_s1)
            arguments
                this mladni.AdniDemographics
                acqdate datetime
                amyloid_s1 table % derived from this.table_amyloid, for single subject
            end

            tf_av45 = strcmp(amyloid_s1.TRACER, 'av45') & amyloid_s1.SUMMARYSUVR_WHOLECEREBNORM >= 1.11;
            tf_fbb = strcmp(amyloid_s1.TRACER, 'fbb') & amyloid_s1.SUMMARYSUVR_WHOLECEREBNORM >= 1.08;
            amyloid_s1.AmyloidStatus = double(tf_av45 | tf_fbb);
            as = this.AmyloidStatus(acqdate, amyloid_s1);
        end
        function as = AmyloidStatusLong(this, acqdate, amyloid_s1)
            arguments
                this mladni.AdniDemographics
                acqdate datetime
                amyloid_s1 table % derived from this.table_amyloid, for single subject
            end

            tf_av45 = strcmp(amyloid_s1.TRACER, 'av45') & amyloid_s1.SUMMARYSUVR_COMPOSITE_REFNORM >= 0.78;
            tf_fbb = strcmp(amyloid_s1.TRACER, 'fbb') & amyloid_s1.SUMMARYSUVR_COMPOSITE_REFNORM >= 0.74;
            amyloid_s1.AmyloidStatus = double(tf_av45 | tf_fbb);
            as = this.AmyloidStatus(acqdate, amyloid_s1);
        end
        function tf = has_fdgfn(this)
            tf = cellfun(@(x) ~isempty(x), this.table_fdg1.FdgFilename);
        end
        function tf = has_ponsvermis(this)
            tf = ~isnan(this.table_fdg1.PonsVermis);
        end
        function tf = has_t1wfn(this)
            tf = cellfun(@(x) ~isempty(x), this.table_fdg1.T1wFilename);
        end
        function [t,t_av45,t_fbb] = table_amyloid(this, varargin)
            [t,t_av45,t_fbb] = this.merge_.table_amyloid(varargin{:});
        end
        function t = table_apoe(this, varargin)
            t = this.merge_.table_apoe(varargin{:});            
        end
        function t = table_cdr(this, varargin)
            t = this.merge_.table_cdr(varargin{:});
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
            t = this.table_paren(t, varargin{:});
        end        
        function t = table_dCDRSB(this, varargin)
            if ~isempty(this.dCDRSB_)
                t = this.dCDRSB_;
                return
            end

            % lazy init
            s = this.subjects_(1);
            t_ = this.merge_;
            t__ = t_(t_.RID == s & ~isnan(t_.FDG), :); % pick subject subtable
            d = max(t__.CDRSB) - min(t__.CDRSB_bl); % scalar dCDRSB
            t = table(s, d, 'VariableNames', {'RID', 'dCDRSB'});
            for ui = 2:length(this.subjects_)
                s = this.subjects_(ui);
                u__ = t_(t_.RID == s & ~isnan(t_.FDG), :); 
                d = max(u__.CDRSB) - min(u__.CDRSB_bl);
                try
                    t = [t; table(s, d, 'VariableNames', {'RID', 'dCDRSB'})]; %#ok<AGROW> 
                catch
                    %fprintf("no reasonable subtable for RID ~ %g\n", this.subjects_(ui))
                end
            end
            this.dCDRSB_ = t;
            t = this.table_paren(t, varargin{:});
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
            t = this.table_paren(t, varargin{:});
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
            t = this.table_paren(t, varargin{:});
        end
        function t = table_fdg(this, varargin)
            %% Table of FDG metrics provided by LONI's downloading services.

            if isempty(this.fdg_)
                this.fdg_ = readtable(this.fdg_proc1_file);
            end
            t = this.fdg_;
            t = this.table_paren(t, varargin{:});
        end
        function t = table_fdg_select_description(this, varargin)
            t = this.table_select_description(this.table_fdg, varargin{:});
        end  
        function t = table_fdg1(this, varargin)
            %%  Build a single table containing mutually consistent data elements for purposes of analyses
            %   of FDG from ADNI.  
            %   - Prefer using csv tables obtained directly from LONI, ensuring csv dates are sufficiently spanning.
            %     See this.table_fdg_proc1_file, generated by LONI at download of FDG imaging, and 
            %     mladni.AdniMerge.*_file for native csv files names.  Collect data acquisition dates specific for each
            %     native csv table, e.g., EXAMDATE, APTESTDT, renaming for disambiguation.  N.B.:
            %     ADNIMERGE*.csv, CDR*.csv, UCBERKELEYFDG*.csv, UCBERKELEYAV45*.csv, UCBERKELEYFBB*.csv, APOERES*.csv,
            %     which are needed for essential covariates age, sex, CDR, FDG references, AV45, FBB and ApoE.
            %   - For data from csv file ADNIMERGE*.csv, replace all EXAMDATES with values from REGISTRY*.csv,
            %     using RID and VISCODE as selectors.  Appears unnecessary since ADNIMERGE_14May2023.csv.
            %   - Select an FDG preprocessing directive, e.g., 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution'
            %   - For each subject, selected by RID:
            %       Use FDG scan acquisition date from this.table_fdg_proc1_file:AcqDate.  
            %       For each scan acquisition date:
            %         For each native csv table: 
            %           manage multiple data collection dates, e.g., for CDR*.csv, USERDATE, USERDATE2, CDDATE, EXAMDATE.
            %           if table is for amyloid:
            %             See also: mladni.AdniDemographics.AmyloidStatus*
            %             (https://github.com/sotiraslab/mladni/blob/main/src/%2Bmladni/AdniDemographics.m function AmyloidStatus*.)

            % cached in memory
            if ~isempty(this.fdg1_)
                fprintf('%s: using cached in memory\n', stackstr)
                t = this.fdg1_;
                if ~isempty(varargin)
                    t = t(varargin{:});
                end       
                return
            end

            % cached on filesystem
            if isfile(this.fdg1_file) && this.reuse_cache
                fprintf('%s: using cached on filesystem\n', stackstr)
                ld = load(this.fdg1_file);
                this.fdg1_ = ld.table_fdg1;
                t = this.fdg1_;
                if ~isempty(varargin)
                    t = t(varargin{:});
                end      
                return
            end
            
            % build caches
            fprintf('%s: building cache\n', stackstr)
            this.fdg1_ = this.buildTableFdg1();
            table_fdg1 = this.fdg1_;
            save(this.fdg1_file, 'table_fdg1');
            t = this.fdg1_;   
            t = this.table_paren(t, varargin{:});         
        end
        function t = table_fdg2(this, varargin)
            %%  this.table_fdg1(this.has_t1wfn, :)

            if ~isempty(this.fdg2_)
                t = this.table_paren(this.fdg2_, varargin{:});
                return
            end
            fdg1 = this.table_fdg1();
            this.fdg2_ = fdg1(this.has_t1wfn, :);
            t = this.table_paren(this.fdg2_, varargin{:});  
        end
        function t = table_fdg3(this, varargin)
            %%  this.table_fdg1(this.has_t1wfn & this.has_ponsvermis, :)
            
            if ~isempty(this.fdg3_)
                t = this.table_paren(this.fdg3_, varargin{:});
                return
            end
            fdg1 = this.table_fdg1();
            this.fdg3_ = fdg1(this.has_t1wfn & this.has_ponsvermis, :);
            t = this.table_paren(this.fdg3_, varargin{:});  
        end
        function t = table_firstscan(this, t_in, varargin)
            %% do not use lazy init with table_fdg()

            t_ = t_in;
            subjects__ = unique(t_.Subject);
            t = t_(strcmp(t_.Subject, subjects__{1}), :); 
            t = t(t.AcqDate == min(t.AcqDate), :); 
            for ti = 2:length(subjects__)
                u = t_(strcmp(t_.Subject, subjects__{ti}), :); % pick subject
                u = u(u.AcqDate == min(u.AcqDate), :); % pick subject's first scan                
                t = [t; u]; %#ok<AGROW> % append 
            end
            this.firstscan_ = t;
            t = this.table_paren(t, varargin{:});
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
            t = this.table_paren(t, varargin{:});
        end     
        function t = table_merge(this, varargin)
            t = this.merge_.table_merge(varargin{:});
        end
        function t = table_neuropath(this, varargin)
            t = this.merge_.table_neuropath(varargin{:});
        end
        function t = table_pet_c3(this, varargin)
            t = this.merge_.table_pet_c3(varargin{:});
        end
        function t = table_pet_qc(this, varargin)
            t = this.merge_.table_pet_qc(varargin{:});
        end
        function t = table_registry(this, varargin)
            t = this.merge_.table_registry(varargin{:});
        end
        function t = table_select_available(~, t_in, varname)
            avail = cellfun(@(x) ~isempty(x), t_in.(varname));
            t = t_in(avail, :);
        end
        function t = table_select_description(this, t_in, varargin)
            t = t_in(strcmp(t_in.Description, this.description), :);
            t = this.table_paren(t, varargin{:});
        end 
        function t = table_strokesum(this, varargin)
            t = this.merge_.table_strokesum(varargin{:});            
        end
        function t = table_tau(this, varargin)
            t = this.merge_.table_tau(varargin{:});
        end
        function t = table_ucberkeleyav1451(this, varargin)
            t = this.merge_.table_ucberkeleyav1451(varargin{:});
        end
        function t = table_ucberkeleyfdg(this, varargin)
            t = this.merge_.table_ucberkeleyfdg(varargin{:});
        end
        function t = table_ucsdvol(this, varargin)
            t = this.merge_.table_ucsdvol(varargin{:});
        end
        function t = table_ucsfsntvol(this, varargin)
            t = this.merge_.table_ucsfsntvol(varargin{:});            
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
        %  subgroup total = 247 + 133 + 148 + 166 + 87 = 781

        function     call(this, cs)
            %% sequentially calls create_subgroups(), create_ic_means()

            arguments
                this mladni.AdniDemographics
                cs logical = false
            end

            error("mladni:NotImplementedError", "AdniDemographics.call:  globs, labels do too much")

            f1 = @this.create_subgroups;
            f2 = @this.create_ic_means;

            aglobs = this.AGLOBS;
            lblg = this.LABELS_G;
            lbl = this.LABELS;

            idx = 1;
            f1(aglobs{idx}, lblg{idx}, cs);
            f2(lbl{idx}, lblg{idx}, cs);
        end
        function     create_ic_means(this, label, labelg, study_design)
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
                study_design {mustBeTextScalar} = 'longitudinal' % cross-sectional
            end
            labelg = strcat('AdniDemographics_create_subgroups_', labelg); % find mat saved by create_subgroups
            
            pwd0 = pushd(this.workdir);
            load(sprintf('%s_%s_fqfns.mat', labelg, study_design)) %#ok<*LOAD> 
            
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
                ic.(this.subgroups{1}).view_qc(icd.(sg{1}));
                ic.(this.subgroups{1}).save_qc(icd.(sg{1}));
            end
            popd(pwd0);
        end 
        function     create_subgroups(this, aglob, labelg, opts)
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
                labelg {mustBeTextScalar} = 'fdg'
                opts.study_design {mustBeTextScalar} = 'longitudinal' % cross-sectional
            end
            
            pwd0 = pushd(this.workdir);
            %to_remove = mybasename(this.fqfns_scans_gt_1y);
            cs = strcmp(opts.study_design, 'cross-sectional');

            for sg = this.subgroups
                table_ = strcat('table_', sg{1});
                TT.(sg{1}) = this.(table_)(cs); % call this.table_*()
                fqfns.(sg{1}) = this.table2fqfns(TT.(sg{1}), globbing=aglob); % call this.table2fqfns()
                
                Filename = fqfns.(sg{1}); % rename for addvars()
                TT.(sg{1}) = addvars(TT.(sg{1}), Filename, 'Before', 1, 'NewVariableNames', {'Filename'});
                %select = cellfun(@(x) ~isempty(x), Filename);
                %disp(TT.(sg{1}))

                writetable(TT.(sg{1}), sprintf('table_%s_%s.csv', sg{1}, opts.study_design), WriteVariableNames=true)
            end

            save(sprintf('AdniDemographics_create_subgroups_%s_%s_fqfns.mat', labelg, opts.study_design), 'fqfns')
            save(sprintf('AdniDemographics_create_subgroups_%s_%s_TT.mat', labelg, opts.study_design), 'TT')
            popd(pwd0);
        end
        function t = table_all(this, cs, T_name)
            %% ...; N = 1544 with late baseline; N = 3478 longitudinal by table_fdg3.

            arguments
                this mladni.AdniDemographics
                cs = false % cross-sectional
                T_name = 'table_fdg3'
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__;
            else
                t = this.table_firstscan(fdg__);
            end
        end
        function t = table_cn(this, cs, T_name, bl_1st)
            %% N = 247 in nifti_files.csv; N = 262 with late baseline; N = 492->479 longitudinal by table_fdg3.
            
            arguments
                this mladni.AdniDemographics
                cs logical = false % cross-sectional
                T_name {mustBeTextScalar} = 'table_fdg3'
                bl_1st logical = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusCS == 0, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusCS == 0, :);
            else
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusCS == 0, :);
                t = this.table_firstscan(t);
            end
        end
        function t = table_preclinical(this, cs, T_name, bl_1st)
            %% N = 133 in nifti_files.csv; N = 140 with late baseline; N = 164->177 longitudinal by table_fdg3.

            arguments
                this mladni.AdniDemographics
                cs = false % cross-sectional
                T_name = 'table_fdg3'
                bl_1st = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusCS == 1, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusCS == 1, :);
            else
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusCS == 1, :);
                t = this.table_firstscan(t);
            end
        end        
        function t = table_cdr_0p5_apos(this, cs, T_name, bl_1st)
            %% N = 385 in nifti_files.csv; N = 460 with late baseline; N = 566->575 longitudinal by table_fdg3.

            arguments
                this mladni.AdniDemographics
                cs = false % cross-sectional
                T_name = 'table_fdg3'
                bl_1st = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL == 0.5 & fdg__.AmyloidStatusCS == 1, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL == 0.5 & fdg__.AmyloidStatusCS == 1, :);
            else
                t = fdg__(fdg__.CDGLOBAL == 0.5 & fdg__.AmyloidStatusCS == 1, :);
                t = this.table_firstscan(t);
            end
        end
        function t = table_cdr_gt_0p5_apos(this, cs, T_name, bl_1st)
            %% N = 87 in nifti_files.csv; N = 149 with late baseline; N = 173->173 longitudinal by table_fdg3.

            arguments
                this mladni.AdniDemographics
                cs = false % cross-sectional
                T_name = 'table_fdg3'
                bl_1st = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL > 0.5 & fdg__.AmyloidStatusCS == 1, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL > 0.5 & fdg__.AmyloidStatusCS == 1, :);
            else
                t = fdg__(fdg__.CDGLOBAL > 0.5 & fdg__.AmyloidStatusCS == 1, :);
                t = this.table_firstscan(t);
            end            
        end
        function t = table_cdr_ge_0p5_aneg(this, cs, T_name, bl_1st)
            %% N = 247; N = 333 with late baseline; N = 564->557 longitudinal by table_fdg3.

            arguments
                this mladni.AdniDemographics
                cs = false % cross-sectional
                T_name = 'table_fdg3'
                bl_1st = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL >= 0.5 & fdg__.AmyloidStatusCS == 0, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL >= 0.5 & fdg__.AmyloidStatusCS == 0, :);
            else
                t = fdg__(fdg__.CDGLOBAL >= 0.5 & fdg__.AmyloidStatusCS == 0, :);
                t = this.table_firstscan(t);
            end            
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
            assert(all(strcmp(last.Subject, first.Subject)), stackstr)
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

            assert(istext(group), stackstr)

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
        merge_
        dCDRSB_
        fdg_
        fdg1_
        fdg2_
        fdg3_
        firstscan_
        fqfns_scans_gt_1y_
        lastscan_
        subjects_
    end

    methods (Access = protected)
        function t = buildTableFdg1(this)
            % Additional sources of FreeSurfer:
            % UCSFFSX.csv: 
            % UCSFFSX51.csv

            % Additional sources of abeta, tau, ptau:
            % UPENNBIOMK.csv:  
            % UPENNBIOMK2.csv
            % UPENNBIOMK3.csv
            % UPENNBIOMK4.csv
            % UPENNBIOMK5.csv
            % UPENNBIOMK6.csv

            t = this.table_fdg(); % enumerates ADNI/bids/derivatives/sub-*/ses-*/*_trc-fdg_*.nii.gz
            t = t(contains(t.Description, this.description), :);
            t.AcqDate = datetime(t.AcqDate);
            sz = size(t.ImageDataID);

            %% fdg filename

            t = addvars(t, cell(sz), Before=1, NewVariableName='T1wFilename');
            t = addvars(t, cell(sz), Before=1, NewVariableName='FdgFilename');

            %% this.table_merge

            t.MergeRid = nan(sz);
            t.MergePtid = cell(sz);
            t.MergeVisCode = cell(sz);
            t.MergeSite = nan(sz);
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
            t.MergeMoca = nan(sz);
            t.MergeImageUid = nan(sz);
            t.MergeVentricles = nan(sz);
            t.MergeHippocampus = nan(sz);
            t.MergeWholeBrain = nan(sz);
            t.MergeIcv = nan(sz);
            t.MergeDx = cell(sz);
            t.MergeMPACCdigit = nan(sz);
            t.MergeMPACCtrailsB = nan(sz);

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
            t.MergeMPACCdigitBl = nan(sz);
            t.MergeMPACCtrailsBBl = nan(sz);
            t.MergeVentriclesBl = nan(sz);
            t.MergeHippocampusBl = nan(sz);
            t.MergeWholeBrainBl = nan(sz);
            t.MergeIcvBl = nan(sz);
            t.MergeMocaBl = nan(sz);
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

            %% this.table_ucberkeleyfdg
            
            t.BerkeleyFdgExamDate = NaT(sz);
            t.Metaroi = nan(sz);
            t.PonsVermis = nan(sz);

            %% this.table_amyloid

            t.AmyloidExamDate = NaT(sz);
            t.TRACER = cell(sz);
            t.AmyloidStatusCS = nan(sz);
            t.AmyloidStatusLong = nan(sz);
            t.SUMMARYSUVR_WHOLECEREBNORM = nan(sz);
            t.SUMMARYSUVR_COMPOSITE_REFNORM = nan(sz);

            %% this.table_tau

            t.TauExamDate = NaT(sz);
            t.META_TEMPORAL_SUVR = nan(sz);
            t.BRAAK1_SUVR = nan(sz);
            t.BRAAK34_SUVR = nan(sz);
            t.BRAAK56_SUVR = nan(sz);
            t.META_TEMPORAL_VOLUME = nan(sz);
            t.BRAAK1_VOLUME = nan(sz);
            t.BRAAK34_VOLUME = nan(sz);
            t.BRAAK56_VOLUME = nan(sz);

            %% this.table_{ucsdvol,ucsfsntvol}

            t.UcsdExamDate = NaT(sz);
            t.UcsdLHippo = nan(sz);
            t.UcsdRHippo = nan(sz);

            t.UcsfExamDate = NaT(sz);
            t.UcsfLHippo = nan(sz);
            t.UcsfRHippo = nan(sz);

            %% this.table_apoe

            t.ApTestDt = NaT(sz);
            t.ApGen1 = nan(sz);
            t.ApGen2 = nan(sz);
            t.ApoE2 = nan(sz);

            %% this.table_neuropath
            
            t.NpBraak = nan(sz);

            %% this.table_cdr

            t.ORIGPROT = cell(sz);
            t.COLPROT = cell(sz);
            t.RID = nan(sz);
            t.SITEID = nan(sz);
            t.VISCODE = cell(sz);
            t.USERDATE = NaT(sz);
            t.USERDATE2 = NaT(sz);
            t.CDCARE = nan(sz);
            t.CDCOMMUN = nan(sz);
            t.CDDATE = NaT(sz);
            t.CDHOME = nan(sz);
            t.CDJUDGE = nan(sz);
            t.CDMEMORY = nan(sz);
            t.CDORIENT = nan(sz);
            t.CDSOURCE = cell(sz);
            t.CDVERSION = cell(sz);
            t.SPID = cell(sz);
            t.CDGLOBAL = nan(sz);
            t.CDEXAMDATE = NaT(sz);
            t.CDRSB = nan(sz);

            %% iterate this.table_fdg.Subject

            %%
            % N of fdg baselines counted by unique RIDs in UCBERKELEYFDG_03_25_22.csv == 1624
            % N of fdg subject folders in ADNI/bids/derivatives == 1659
            % >> size(unique(this.subjects)) % ~ this.table_fdg.Subject ~ \d{3}_S_\d{4}
            %         1660           1
            % >> size(unique(this.table_merge.PTID)) % ~ ADNIMERGE.csv.PTID ~ \d{3}_S_\d{4}
            %         2430           1
            % >> size(unique(this.table_merge.RID)) % ~ ADNIMERGE.csv.PTID ~ \d+
            %         2430           1      
            % >> size(intersect(this.table_merge.PTID, this.table_fdg.Subject))
            %         1656           1            

            N_skipped_T1w = 0;
            N_skipped_merge = 0;
            N_skipped_berkeley_metaroi = 0;
            N_skipped_berkeley_ponsvermis = 0;
            N_skipped_ucsd_hippo = 0;
            N_skipped_ucsf_hippo = 0;
            N_skipped_apoe = 0;
            N_skipped_cdr = 0;
            N_skipped_amyloid = 0;

            %%

            for sub = asrow(this.subjects)
                try
                    t_s1 = t(strcmp(t.Subject, sub{1}), :); % consider:  t.Subject -> this.subjects

                    %% create subtables specific to subject & data source
                    
                    try
                        merge_s1 = this.table_merge(strcmp(this.table_merge.PTID, sub{1}), ':');
                        rid = merge_s1.RID(1);
                    catch 
                        fprintf('%s: table_merge(%s) is missing\n', stackstr, sub{1});
                        % continue
                    end
                    
                    try
                        berkeley_fdg_s1 = this.table_ucberkeleyfdg(this.table_ucberkeleyfdg.RID == rid, ':');                        
                        berkeley_metaroi_s1 = this.table_ucberkeleyfdg(this.table_ucberkeleyfdg.RID == rid, ':');
                        berkeley_metaroi_s1 = berkeley_metaroi_s1(contains(berkeley_metaroi_s1.ROINAME, 'metaroi'), :);    
                        berkeley_ponsvermis_s1 = this.table_ucberkeleyfdg(this.table_ucberkeleyfdg.RID == rid, ':');
                        berkeley_ponsvermis_s1 = berkeley_ponsvermis_s1(contains(berkeley_ponsvermis_s1.ROINAME, 'pons-vermis'), :);
                    catch
                        % fprintf('%s: table_ucberkeleyfdg(%s) is missing\n', stackstr, sub{1});
                        % continue
                    end

                    try
                        t_amy = this.table_amyloid;
                        select = t_amy.RID == rid;
                        assert(islogical(select), stackstr)
                        amyloid_s1 = t_amy(select, :); % Incipient BUG ?
                    catch ME
                        disp(ME)
                        % fprintf('%s: table_amyloid(%s) is missing\n', stackstr, sub{1}); % too frequent to fprint
                        % continue
                    end

                    try
                        tau_s1 = this.table_ucsfsntvol(this.table_ucsfsntvol.RID == rid, ':');
                    catch 
                    end

                    try
                        ucsd_hippo_s1 = this.table_ucsdvol(this.table_ucsdvol.RID == rid, ':');
                        ucsf_hippo_s1 = this.table_ucsfsntvol(this.table_ucsfsntvol.RID == rid, ':');
                    catch 
                    end

                    try
                        apoe_s1 = this.table_apoe(this.table_apoe.RID == rid, ':');
                    catch 
                    end

                    try
                        t_np = this.table_neuropath;
                        select = t_np.RID == rid;
                        assert(islogical(select))
                        neuropath_s1 = t_np(select, :); % Incipient BUG ?
                    catch 
                    end

                    try
                        cdr_s1 = this.table_cdr(this.table_cdr.RID == rid, ':');
                    catch 
                        % fprintf('%s: table_cdr(%s) is missing\n', stackstr, sub{1});
                        % continue
                    end

                    %% iterate single subject's t_s1.AcqDate
    
                    numrows_s1 = size(t_s1, 1);
                    for acq = 1:numrows_s1
                        acqdate = t_s1.AcqDate(acq);
                        acqdate.Format = 'default';

                        try
                            t_s1.FdgFilename{acq} = this.table2fdgfilename(t_s1(acq, :));
                            t_s1.T1wFilename{acq} = this.table2t1wfilename(t_s1(acq, :));
                            
                            sep = this.filenames2sep(t_s1.FdgFilename{acq}, t_s1.T1wFilename{acq});
                            if duration(sep) > duration(this.datetime_separation_tol)
                                throw(MException('mladni:datetime_separation_tol_exceeded', ...
                                    'FdgFilename, T1wFilename datetime separation = %s years', years(sep)))
                            end
                        catch ME
                            N_skipped_T1w = N_skipped_T1w + 1;
                            if ~contains(ME.message, 'table2t1wfilename') && ~contains(ME.message, 'table2fdgfilename')
                                disp(ME.message)
                            end
                        end

                        try
                            [sep,merge_near] = min(abs(acqdate - merge_s1.EXAMDATE));
                            if duration(sep) > duration(this.datetime_separation_tol)
                                N_skipped_merge = N_skipped_merge + 1;
                                throw(MException('mladni:datetime_separation_tol_exceeded', ...
                                    'merge_s1 datetime separation = %s years', years(sep)))
                            end
                            t_s1(acq, 'MergeRid') = merge_s1(merge_near, 'RID');
                            t_s1(acq, 'MergePtid') = merge_s1(merge_near, 'PTID');
                            t_s1(acq, 'MergeVisCode') = merge_s1(merge_near, 'VISCODE');
                            t_s1(acq, 'MergeSite') = merge_s1(merge_near, 'SITE');
                            t_s1(acq, 'MergeExamDate') = merge_s1(merge_near, 'EXAMDATE'); % ------- Merge EXAMDATE -------
                            t_s1(acq, 'MergeDxBl') = merge_s1(merge_near, 'DX_bl');
                            t_s1(acq, 'MergeAge') = merge_s1(merge_near, 'AGE');
                            t_s1(acq, 'MergePtGender') = merge_s1(merge_near, 'PTGENDER');
                            t_s1(acq, 'MergePtEducat') = merge_s1(merge_near, 'PTEDUCAT');
                            t_s1(acq, 'MergePtEthCat') = merge_s1(merge_near, 'PTETHCAT');
                            t_s1(acq, 'MergePtRacCat') = merge_s1(merge_near, 'PTRACCAT');
                            t_s1(acq, 'MergePtMarry') = merge_s1(merge_near, 'PTMARRY');
                            t_s1(acq, 'MergeApoE4') = merge_s1(merge_near, 'APOE4');
                            t_s1(acq, 'MergeFdg') = merge_s1(merge_near, 'FDG');
                            t_s1(acq, 'MergePib') = merge_s1(merge_near, 'PIB');
                            t_s1(acq, 'MergeAv45') = merge_s1(merge_near, 'AV45');
                            t_s1(acq, 'MergeAbeta') = merge_s1(merge_near, 'ABETA');
                            t_s1(acq, 'MergeTau') = merge_s1(merge_near, 'TAU');
                            t_s1(acq, 'MergePTau') = merge_s1(merge_near, 'PTAU');
                            t_s1(acq, 'MergeCdrsb') = merge_s1(merge_near, 'CDRSB');
                            t_s1(acq, 'MergeMmse') = merge_s1(merge_near, 'MMSE');
                            t_s1(acq, 'MergeMoca') = merge_s1(merge_near, 'MOCA');
                            t_s1(acq, 'MergeImageUid') = merge_s1(merge_near, 'IMAGEUID');
                            t_s1(acq, 'MergeVentricles') = merge_s1(merge_near, 'Ventricles');
                            t_s1(acq, 'MergeHippocampus') = merge_s1(merge_near, 'Hippocampus');
                            t_s1(acq, 'MergeWholeBrain') = merge_s1(merge_near, 'WholeBrain');
                            t_s1(acq, 'MergeIcv') = merge_s1(merge_near, 'ICV');
                            t_s1(acq, 'MergeDx') = merge_s1(merge_near, 'DX');
                            t_s1(acq, 'MergeMPACCdigit') = merge_s1(merge_near, 'mPACCdigit');
                            t_s1(acq, 'MergeMPACCtrailsB') = merge_s1(merge_near, 'mPACCtrailsB');
                            t_s1(acq, 'MergeExamDateBl') = merge_s1(merge_near, 'EXAMDATE_bl');
                            t_s1(acq, 'MergeCdrsbBl') = merge_s1(merge_near, 'CDRSB_bl');
                            t_s1(acq, 'MergeAdas11Bl') = merge_s1(merge_near, 'ADAS11_bl');
                            t_s1(acq, 'MergeAdas13Bl') = merge_s1(merge_near, 'ADAS13_bl');
                            t_s1(acq, 'MergeAdasq4Bl') = merge_s1(merge_near, 'ADASQ4_bl');
                            t_s1(acq, 'MergeMmseBl') = merge_s1(merge_near, 'MMSE_bl');
                            t_s1(acq, 'MergeRavltImmediateBl') = merge_s1(merge_near, 'RAVLT_immediate_bl');
                            t_s1(acq, 'MergeRavltLearningBl') = merge_s1(merge_near, 'RAVLT_learning_bl');
                            t_s1(acq, 'MergeRavltForgettingBl') = merge_s1(merge_near, 'RAVLT_forgetting_bl');
                            t_s1(acq, 'MergeRavltPercForgettingBl') = merge_s1(merge_near, 'RAVLT_perc_forgetting_bl');
                            t_s1(acq, 'MergeMPACCdigitBl') = merge_s1(merge_near, 'mPACCdigit_bl');
                            t_s1(acq, 'MergeMPACCtrailsBBl') = merge_s1(merge_near, 'mPACCtrailsB_bl');
                            t_s1(acq, 'MergeVentriclesBl') = merge_s1(merge_near, 'Ventricles_bl');
                            t_s1(acq, 'MergeHippocampusBl') = merge_s1(merge_near, 'Hippocampus_bl');
                            t_s1(acq, 'MergeWholeBrainBl') = merge_s1(merge_near, 'WholeBrain_bl');
                            t_s1(acq, 'MergeIcvBl') = merge_s1(merge_near, 'ICV_bl');
                            t_s1(acq, 'MergeMocaBl') = merge_s1(merge_near, 'MOCA_bl');
                            t_s1(acq, 'MergeEcogPtMemBl') = merge_s1(merge_near, 'EcogPtMem_bl');
                            t_s1(acq, 'MergeEcogPtLangBl') = merge_s1(merge_near, 'EcogPtLang_bl');
                            t_s1(acq, 'MergeEcogPtVisspatBl') = merge_s1(merge_near, 'EcogPtVisspat_bl');
                            t_s1(acq, 'MergeEcogPtPlanBl') = merge_s1(merge_near, 'EcogPtPlan_bl');
                            t_s1(acq, 'MergeEcogPtOrganBl') = merge_s1(merge_near, 'EcogPtOrgan_bl');
                            t_s1(acq, 'MergeEcogPtDivattBl') = merge_s1(merge_near, 'EcogPtDivatt_bl');
                            t_s1(acq, 'MergeEcogPtTotalBl') = merge_s1(merge_near, 'EcogPtTotal_bl');
                            t_s1(acq, 'MergeEcogSPMemBl') = merge_s1(merge_near, 'EcogSPMem_bl');
                            t_s1(acq, 'MergeEcogSPLangBl') = merge_s1(merge_near, 'EcogSPLang_bl');
                            t_s1(acq, 'MergeEcogSPVisspatBl') = merge_s1(merge_near, 'EcogSPVisspat_bl');
                            t_s1(acq, 'MergeEcogSPPlanBl') = merge_s1(merge_near, 'EcogSPPlan_bl');
                            t_s1(acq, 'MergeEcogSPOrganBl') = merge_s1(merge_near, 'EcogSPOrgan_bl');
                            t_s1(acq, 'MergeEcogSPDivattBl') = merge_s1(merge_near, 'EcogSPDivatt_bl');
                            t_s1(acq, 'MergeEcogSPTotalBl') = merge_s1(merge_near, 'EcogSPTotal_bl');
                            t_s1(acq, 'MergeAbetaBl') = merge_s1(merge_near, 'ABETA_bl');
                            t_s1(acq, 'MergeTauBl') = merge_s1(merge_near, 'TAU_bl');
                            t_s1(acq, 'MergePTauBl') = merge_s1(merge_near, 'PTAU_bl');
                            t_s1(acq, 'MergeFdgBl') = merge_s1(merge_near, 'FDG_bl');
                            t_s1(acq, 'MergePibBl') = merge_s1(merge_near, 'PIB_bl');
                            t_s1(acq, 'MergeAv45Bl') = merge_s1(merge_near, 'AV45_bl');
                            t_s1(acq, 'MergeYearsBl') = merge_s1(merge_near, 'Years_bl');
                            t_s1(acq, 'MergeMonthBl') = merge_s1(merge_near, 'Month_bl');
                        catch ME
                            if contains(ME.message, 'datetime separation')
                                disp(ME.message)
                            end
                        end

                        try
                            [sep,berkeley_near] = min(abs(acqdate - berkeley_metaroi_s1.EXAMDATE));
                            if duration(sep) > duration(this.datetime_separation_tol)
                                N_skipped_berkeley_metaroi = N_skipped_berkeley_metaroi + 1;                                
                                throw(MException('mladni:datetime_separation_tol_exceeded', ...
                                    'berkeley_metaroi_s1 datetime separation = %s years', years(sep)))
                            end
                            t_s1(acq, 'Metaroi') = berkeley_metaroi_s1(berkeley_near, 'MEAN');   

                            [sep,berkeley_near] = min(abs(acqdate - berkeley_ponsvermis_s1.EXAMDATE));
                            if duration(sep) > duration(this.datetime_separation_tol)
                                N_skipped_berkeley_ponsvermis = N_skipped_berkeley_ponsvermis + 1;
                                throw(MException('mladni:datetime_separation_tol_exceeded', ...
                                    'berkeley_ponsvermis_s1 datetime separation = %s years', years(sep)))
                            end
                            t_s1(acq, 'PonsVermis') = berkeley_ponsvermis_s1(berkeley_near, 'MEAN');

                            t_s1(acq, 'BerkeleyFdgExamDate') = berkeley_fdg_s1(berkeley_near, 'EXAMDATE');
                        catch ME
                            if contains(ME.message, 'datetime separation')
                                disp(ME.message)
                            end
                        end

                        try
                            assert(~isempty(amyloid_s1), ...
                                '%s: amyloid_s1 empty for sub %s acqdate %s', stackstr, sub{1}, acqdate)
                            [~,amyloid_near] = min(abs(acqdate - amyloid_s1.EXAMDATE));
                            t_s1(acq, 'AmyloidExamDate') = amyloid_s1(amyloid_near, 'EXAMDATE');
                            t_s1(acq, 'TRACER') = amyloid_s1(amyloid_near, 'TRACER');
                            t_s1(acq, 'SUMMARYSUVR_WHOLECEREBNORM') = amyloid_s1(amyloid_near, 'SUMMARYSUVR_WHOLECEREBNORM');
                            t_s1(acq, 'SUMMARYSUVR_COMPOSITE_REFNORM') = amyloid_s1(amyloid_near, 'SUMMARYSUVR_COMPOSITE_REFNORM');
                            t_s1{acq, 'AmyloidStatusCS'} = this.AmyloidStatusCS(acqdate, amyloid_s1);
                            t_s1{acq, 'AmyloidStatusLong'} = this.AmyloidStatusLong(acqdate, amyloid_s1);
                        catch ME
                            N_skipped_amyloid = N_skipped_amyloid + 1;
                            if ~contains(ME.message, 'amyloid_s1 empty for sub')
                                disp(ME.message)
                            end
                        end

                        try
                            [sep,tau_near] = min(abs(acqdate - tau_s1.EXAMDATE));
                            if duration(sep) > duration(this.datetime_separation_tol)
                                throw(MException('mladni:datetime_separation_tol_exceeded', ...
                                    'tau_s1 datetime separation = %s years', years(sep)))
                            end
                            t_s1(acq, 'TauExamDate') = tau_s1(tau_near, 'EXAMDATE');
                            t_s1(acq, 'META_TEMPORAL_SUVR') = tau_s1(tau_near, 'META_TEMPORAL_SUVR');
                            t_s1(acq, 'BRAAK1_SUVR') = tau_s1(tau_near, 'BRAAK1_SUVR');
                            t_s1(acq, 'BRAAK34_SUVR') = tau_s1(tau_near, 'BRAAK34_SUVR');
                            t_s1(acq, 'BRAAK56_SUVR') = tau_s1(tau_near, 'BRAAK56_SUVR');
                            t_s1(acq, 'META_TEMPORAL_VOLUME') = tau_s1(tau_near, 'META_TEMPORAL_VOLUME');
                            t_s1(acq, 'BRAAK1_VOLUME') = tau_s1(tau_near, 'BRAAK1_VOLUME');
                            t_s1(acq, 'BRAAK34_VOLUME') = tau_s1(tau_near, 'BRAAK34_VOLUME');
                            t_s1(acq, 'BRAAK56_VOLUME') = tau_s1(tau_near, 'BRAAK56_VOLUME');
                        catch 
                            % BRAAK staging is post-mortem, so tau_s1 is mostly empty
                        end

                        try
                            [sep,ucsd_near] = min(abs(acqdate - ucsd_hippo_s1.EXAMDATE));
                            if duration(sep) > duration(this.datetime_separation_tol)
                                N_skipped_ucsd_hippo = N_skipped_ucsd_hippo + 1;
                                throw(MException('mladni:datetime_separation_tol_exceeded', ...
                                    'ucsd_hippo_s1 datetime separation = %s years', years(sep)))
                            end
                            t_s1(acq, 'UcsdExamDate') = ucsd_hippo_s1(ucsd_near, 'EXAMDATE');
                            t_s1(acq, 'UcsdLHippo') = ucsd_hippo_s1(ucsd_near, 'LHIPPOC');
                            t_s1(acq, 'UcsdRHippo') = ucsd_hippo_s1(ucsd_near, 'RHIPPOC');
                        catch ME
                            %if contains(ME.message, 'datetime separation')
                            %    disp(ME.message)
                            %end
                        end

                        try
                            [sep,ucsf_near] = min(abs(acqdate - ucsf_hippo_s1.EXAMDATE));
                            if duration(sep) > duration(this.datetime_separation_tol)
                                N_skipped_ucsf_hippo = N_skipped_ucsf_hippo + 1;
                                throw(MException('mladni:datetime_separation_tol_exceeded', ...
                                    'ucsf_hippo_s1 datetime separation = %s years', years(sep)))
                            end
                            t_s1(acq, 'UcsfExamDate') = ucsf_hippo_s1(ucsf_near, 'EXAMDATE');
                            t_s1(acq, 'UcsfLHippo') = ucsf_hippo_s1(ucsf_near, 'LEFTHIPPO');
                            t_s1(acq, 'UcsfRHippo') = ucsf_hippo_s1(ucsf_near, 'RIGHTHIPPO');
                        catch ME
                            %if contains(ME.message, 'datetime separation')
                            %    disp(ME.message)
                            %end
                        end

                        try
                            [sep,apoe_near] = min(abs(acqdate - apoe_s1.APTESTDT));
                            if duration(sep) > duration(this.datetime_separation_tol)
                                N_skipped_apoe = N_skipped_apoe + 1;
                                throw(MException('mladni:datetime_separation_tol_exceeded', ...
                                    'apoe_s1 datetime separation = %s years', years(sep)))
                            end
                            t_s1(acq, 'ApTestDt') = apoe_s1(apoe_near, 'APTESTDT');
                            t_s1(acq, 'ApGen1') = apoe_s1(apoe_near, 'APGEN1');
                            t_s1(acq, 'ApGen2') = apoe_s1(apoe_near, 'APGEN2');
                            ApoE2 = (apoe_s1{apoe_near, 'APGEN1'} == 2) + ...
                                (apoe_s1{apoe_near, 'APGEN2'} == 2);
                            t_s1(acq, 'ApoE2') = table(ApoE2);
                        catch ME
                            %if contains(ME.message, 'datetime separation')
                            %    disp(ME.message)
                            %end
                        end

                        try
                            t_s1(acq, 'NpBraak') = neuropath_s1(1, 'NPBRAAK');
                        catch 
                            % BRAAK staging is post-mortem, so neuropath_s1 is mostly empty
                        end

                        try
                            % impute nans
                            nats_ = isnat(cdr_s1.EXAMDATE);
                            nats_u_ = isnat(cdr_s1.USERDATE);
                            cdr_s1.EXAMDATE(nats_ & ~nats_u_) = cdr_s1.USERDATE(nats_ & ~nats_u_);
                            % iterate
                            nats_ = isnat(cdr_s1.EXAMDATE);
                            nats_u_ = isnat(cdr_s1.USERDATE2);
                            cdr_s1.EXAMDATE(nats_ & ~nats_u_) = cdr_s1.USERDATE2(nats_ & ~nats_u_);
                            % iterate
                            nats_ = isnat(cdr_s1.EXAMDATE);
                            nats_date_ = isnat(cdr_s1.CDDATE);
                            cdr_s1.EXAMDATE(nats_ & ~nats_date_) = cdr_s1.CDDATE(nats_ & ~nats_date_);

                            [sep,cdr_near] = min(abs(acqdate - cdr_s1.EXAMDATE));
                            if duration(sep) > duration(this.datetime_separation_tol)
                                N_skipped_cdr = N_skipped_cdr + 1;
                                throw(MException('mladni:datetime_separation_tol_exceeded', ...
                                    'cdr_s1 datetime separation = %s years', years(sep)))
                            end
                            t_s1(acq, 'ORIGPROT') = cdr_s1(cdr_near, 'ORIGPROT');
                            t_s1(acq, 'COLPROT') = cdr_s1(cdr_near, 'COLPROT');
                            t_s1(acq, 'RID') = cdr_s1(cdr_near, 'RID');
                            t_s1(acq, 'SITEID') = cdr_s1(cdr_near, 'SITEID'); 
                            t_s1(acq, 'VISCODE') = cdr_s1{cdr_near, 'VISCODE'};
                            t_s1{acq, 'USERDATE'} = cdr_s1{cdr_near, 'USERDATE'};
                            t_s1{acq, 'USERDATE2'} = cdr_s1{cdr_near, 'USERDATE2'};
                            t_s1(acq, 'CDCARE') = cdr_s1(cdr_near, 'CDCARE');
                            t_s1(acq, 'CDCOMMUN') = cdr_s1(cdr_near, 'CDCOMMUN');
                            t_s1(acq, 'CDDATE') = cdr_s1(cdr_near, 'CDDATE');
                            t_s1(acq, 'CDHOME') = cdr_s1(cdr_near, 'CDHOME');
                            t_s1(acq, 'CDJUDGE') = cdr_s1(cdr_near, 'CDJUDGE');
                            t_s1(acq, 'CDMEMORY') = cdr_s1(cdr_near, 'CDMEMORY');
                            t_s1(acq, 'CDORIENT') = cdr_s1(cdr_near, 'CDORIENT');
                            t_s1(acq, 'CDSOURCE') = cdr_s1(cdr_near, 'CDSOURCE');
                            t_s1(acq, 'CDVERSION') = cdr_s1(cdr_near, 'CDVERSION');
                            t_s1(acq, 'SPID') = cdr_s1(cdr_near, 'SPID');
                            t_s1(acq, 'CDGLOBAL') = cdr_s1(cdr_near, 'CDGLOBAL');  
                            t_s1{acq, 'CDEXAMDATE'} = cdr_s1{cdr_near, 'EXAMDATE'}; % ------- CDR EXAMDATE -------
                            t_s1(acq, 'CDRSB') = cdr_s1(cdr_near, 'CDRSB'); 
                        catch ME
                            if contains(ME.message, 'datetime separation')
                                disp(ME.message)
                            end
                        end
                    end
                    
                    t(strcmp(t.Subject, sub{1}), :) = t_s1;
                catch ME
                    handwarning(ME)
                end
            end

            fprintf('%s: N_skipped_T1w = %g\n', stackstr, N_skipped_T1w)
            fprintf('%s: N_skipped_merge = %g\n', stackstr, N_skipped_merge)
            fprintf('%s: N_skipped_berkeley_metaroi = %g\n', stackstr, N_skipped_berkeley_metaroi)
            fprintf('%s: N_skipped_berkeley_ponsvermis = %g\n', stackstr, N_skipped_berkeley_ponsvermis)
            fprintf('%s: N_skipped_ucsd_hippo = %g\n', stackstr, N_skipped_ucsd_hippo)
            fprintf('%s: N_skipped_ucsf_hippo = %g\n', stackstr, N_skipped_ucsf_hippo)
            fprintf('%s: N_skipped_apoe = %g\n', stackstr, N_skipped_apoe)
            fprintf('%s: N_skipped_cdr = %g\n', stackstr, N_skipped_cdr)
            fprintf('%s: N_skipped_amyloid = %g\n', stackstr, N_skipped_amyloid)
        end
    end

    methods % (Access = protected)
        function [fqfn,acqdat] = table2fdgfilename(this, T)
            %% T1w NIfTI must have been linearly registrable to T1w.  Does not check for ponsvermis availability.

            assert(size(T,1) == 1, strcat('assert fail:', stackstr))
            descr = this.description_acro;

            % sub
            re = regexp(T.Subject{1}, '(?<siteid>\d{3})_S_(?<rid>\d{4})', 'names');
            sub = sprintf('sub-%sS%s', re.siteid, re.rid);

            % candidate filenames
            fdgfn = sprintf('*_trc-FDG_proc-%s_orient-rpi_pet.nii.gz', descr);
            globbed = glob(fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', sub, 'ses-*', 'pet', fdgfn));
            assert(~isempty(globbed), strcat(stackstr, sprintf(': %s not found', fdgfn)))
            acqdate = NaT(size(globbed));
            for ig = 1:length(globbed)
                re1 = regexp(globbed{ig}, ...
                    sprintf("\\S+/sub-\\d{3}S\\d{4}_ses-(?<dt>\\d{14})_trc-FDG_proc-%s_orient-rpi_pet.nii.gz", descr), ...
                    'names');
                acqdate(ig) = datetime(re1.dt, InputFormat='yyyyMMddHHmmss');
            end
            T1 = table(globbed, acqdate);
            [~,select] = min(abs(T1.acqdate - T.AcqDate));
            fqfn = T1.globbed{select};
            acqdat = T1.acqdate(select);
            assert(abs(duration(acqdat - T.AcqDate)) < this.datetime_separation_tol, strcat('assert fail:', stackstr))
        end
        function [fqfn,acqdat] = table2t1wfilename(this, T)
            %% T1w must have been warpable to MNI152

            assert(size(T,1) == 1, strcat('assert fail:', stackstr))

            % sub
            re = regexp(T.Subject{1}, '(?<siteid>\d{3})_S_(?<rid>\d{4})', 'names');
            sub = sprintf('sub-%sS%s', re.siteid, re.rid);

            % candidate filenames
            t1wfqfn = fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', sub, 'ses-*', 'pet', ...
                'sub-*_ses-*_acq-*_proc-*_orient-rpi_T1w_brain_2Warp.nii.gz');
            globbed = glob(t1wfqfn);
            assert(~isempty(globbed), strcat(stackstr, sprintf(': %s not found', t1wfqfn)))
            acqdate = NaT(size(globbed));
            for ig = 1:length(globbed)
                re1 = regexp(globbed{ig}, ...
                    '\S+/sub-\d{3}S\d{4}_ses-(?<dt>\d{14})_\S+orient-rpi_T1w_brain_2Warp.nii.gz', ...
                    'names');
                acqdate(ig) = datetime(re1.dt, InputFormat='yyyyMMddHHmmss');
            end
            T1 = table(globbed, acqdate);
            [~,select] = min(abs(T1.acqdate - T.AcqDate));
            fqfn = T1.globbed{select};
            acqdat = T1.acqdate(select);
            assert(abs(duration(acqdat - T.AcqDate)) < this.datetime_separation_tol, strcat('assert fail:', stackstr))
        end
        function sep = filenames2sep(~, fn, fn1)
            re = regexp(fn, '\S+sub-\d{3}S\d{4}_ses-(?<dt>\d{14})_\S+.nii.gz', 'names');
            dt = datetime(re.dt, InputFormat='yyyyMMddHHmmss');
            re1 = regexp(fn1, '\S+sub-\d{3}S\d{4}_ses-(?<dt>\d{14})_\S+.nii.gz', 'names');
            dt1 = datetime(re1.dt, InputFormat='yyyyMMddHHmmss');
            sep = abs(duration(dt - dt1));
        end
    end

    methods (Static, Access = protected)
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
