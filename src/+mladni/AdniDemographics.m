classdef AdniDemographics < handle
    %% queries "all_FDG_20211123_12_14_2021.csv"
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
    end

    properties
        description = 'Coreg, Avg, Standardized Image and Voxel Size'
    end

    properties (Constant)
        categories = {'CN', 'EMCI', 'MCI', 'LMCI', 'SMC', 'AD'};
        viscode2_months = [0 6 12 18 24 36 48 54 60 66 72 78 84 90 96 108 120 126 132 138 144 150];
    end

    properties (Dependent)
        cdr_file
        fdg_orig_file
        fdg_proc1_file
        fdg_proc2_file
        fdg1_file
        fdgproc_filenames_file
        merge_file
        mpr_meta_file
        mri_quality_file
        mri_quality_adni3_file
        pet_c3_file
        pet_meta_adni1_file
        pet_meta_adnigo2_file
        pet_meta_adni3_file
        pet_meta_list_file
        pet_qc_file
        subjects
        t1_filenames_file
    end

    methods

        %% GET

        function g = get.cdr_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "CDR.csv");
            % unique RID ~ 3418
        end
        function g = get.fdg_orig_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "staging", "FDG_original_20211118", "all_FDG_20211118_11_18_2021.csv");
            % unique Subject ~ 1662
        end
        function g = get.fdg_proc1_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "staging", "FDG_processed_20211123", "all_FDG_20211123_12_14_2021.csv");
            % unique Subject ~ 1660, unique ImageDataID ~ 14358
        end
        function g = get.fdg_proc2_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "staging", "FDG_processed_20220104", "prev_failed_FDG_20220104_1_04_2022.csv");
            % unique Subject ~ 5, included in fdg_proc1_file
        end
        function g = get.fdg1_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "staging", "FDG_processed_20220104", "mladni_AdniDemographics_fdg1.csv");
        end
        function g = get.fdgproc_filenames_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "bids", "rawdata", "rosters", "fdgproc_filenames.csv");
        end
        function g = get.merge_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "ADNIMERGE.csv");
            % unique PTID ~ 1855, unique IMAGEUID ~ 13923
        end
        function g = get.mpr_meta_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "MPRAGEMETA.csv");
            % unique SubjectID ~ 2519, unique ImageUID ~ 37574
        end
        function g = get.mri_quality_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "MRIQUALITY.csv");
            % unique RID ~ 840, unique LONIUID ~ 11573
        end
        function g = get.mri_quality_adni3_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "MAYOADIRL_MRI_QUALITY_ADNI3.csv");
            % unique PTID ~ 1035, unique LONI_IMAGE ~ 21744
        end
        function g = get.pet_c3_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "PETC3.csv");
            % unique RID ~ 395, unique LONIUID ~ 407
        end
        function g = get.pet_meta_adni1_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "PETMETA_ADNI1.csv");
            % unique RID ~ 420
        end
        function g = get.pet_meta_adnigo2_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "PETMETA_ADNIGO2.csv");
            % unique RID ~ 1212
        end
        function g = get.pet_meta_adni3_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "PETMETA3.csv");
            % unique RID ~ 641
        end
        function g = get.pet_meta_list_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "PET_META_LIST.csv");
            % unique Subject ~ 2028, unique ImageID ~ 50843
        end
        function g = get.pet_qc_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "studydata", "PETQC.csv");
            % unique RID ~ 1413, unique LONIUID ~ 3950
        end
        function g = get.subjects(this)
            if isempty(this.subjects_)
                tbl_fdg = this.table_fdg();
                this.subjects_ = unique(tbl_fdg.Subject);
            end
            g = this.subjects_;
        end
        function g = get.t1_filenames_file(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "bids", "rawdata", "rosters", "t1_filenames.csv");
        end

        %%

        function this = AdniDemographics(varargin)
            warning('off', 'MATLAB:table:ModifiedAndSavedVarnames') 
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
            if isempty(this.fdg1_)
                this.fdg1_ = this.buildTableFdg1();
            end
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
            % lazy init
            if isempty(this.firstscan_)
                t_ = this.table_fdg_select_description();
                t = t_(strcmp(t_.Subject, this.subjects_{1}), :); 
                t = t(t.AcqDate == min(t.AcqDate), :); 
                for ti = 2:length(this.subjects_)
                    u = t_(strcmp(t_.Subject, this.subjects_{ti}), :); % pick subject
                    v = u(u.AcqDate == min(u.AcqDate), :); % pick first scan                
                    t = [t; v]; %#ok<AGROW> % append 
                end
                this.firstscan_ = t;
            end
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_lastscan(this, varargin)
            % lazy init
            if isempty(this.lastscan_)
                t_ = this.table_fdg_select_description();
                t = t_(strcmp(t_.Subject, this.subjects_{1}), :); 
                t = t(t.AcqDate == max(t.AcqDate), :); 
                for ti = 2:length(this.subjects_)
                    u = t_(strcmp(t_.Subject, this.subjects_{ti}), :); % pick subject
                    v = u(u.AcqDate == max(u.AcqDate), :); % pick first scan                
                    t = [t; v]; %#ok<AGROW> % append 
                end
                this.lastscan_ = t;
            end
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_merge(this, varargin)
            if isempty(this.merge_)
                this.merge_ = readtable(this.merge_file);
            end
            t = this.merge_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_mpr_meta(this, varargin)
            if isempty(this.mpr_meta_)
                this.mpr_meta_ = readtable(this.mpr_meta_file);
            end
            t = this.mpr_meta_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
            
            %t = readtable(this.mpr_meta_file);
            %select_subs = contains(t.SubjectID, this.subjects);
            %t = t(select_subs, :);   

            %select_scans = contains(t.Type, 'Processed');
            %t = t(select_scans, :);

            %select_scans = contains(t.Description, 'correct', 'IgnoreCase', true);
            %t = t(select_scans, :);

            %select_scans = contains(t.Sequence, 'mpr', 'IgnoreCase', true) | ...
            %               contains(t.Sequence, 'rage', 'IgnoreCase', true) | ...
            %               contains(t.Sequence, 'mt1', 'IgnoreCase', true) | ...
            %               contains(t.Sequence, 'ir-fspgr', 'IgnoreCase', true);
            %t = t(select_scans, :);

            %unwanted = ...
            %    {'localizer' 'calibration' 'loc' 'mapping' 'surv' 'fgre' 'scout' 'smartbrain' 't2' 'fmri' ...
            %     'calibration'};
            %select_scans = ~contains(t.Sequence, unwanted, 'IgnoreCase', true);
            %t = t(select_scans, :);

            %this.mpr_meta_ = t;
        end
        function t = table_mri_quality(this, varargin)
            if isempty(this.mri_quality_)
                this.mri_quality_ = readtable(this.mri_quality_file);
            end
            t = this.mri_quality_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function t = table_mri_quality_adni3(this, varargin)
            if isempty(this.mri_quality_adni3_)
                this.mri_quality_adni3_ = readtable(this.mri_quality_adni3_file);
            end
            t = this.mri_quality_adni3_;
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
        function t = table_t1_filenames(this, varargin)
            if isempty(this.t1_filenames_)
                this.t1_filenames_ = readtable(this.t1_filenames_file);
            end
            t = this.t1_filenames_;
            if ~isempty(varargin)
                t = t(varargin{:});
            end
        end
        function write_imageids(~, tbl, fname)
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
            last = this.table_lastscan();
            first = this.table_firstscan();
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
        cdr_
        fdg_
        fdg1_
        fdgproc_filenames_
        firstscan_
        lastscan_
        merge_
        mpr_meta_
        mri_quality_
        mri_quality_adni3_
        pet_c3_
        pet_qc_
        subjects_
        t1_filenames_
    end

    methods (Access = protected)
        function t = buildTableFdg1(this)
            
            t = this.table_fdg();
            sz = size(t.AcqDate);

            t.AcqDate = datetime(t.AcqDate);
            t.MergeVisCode = cell(sz);
            t.MergeExamDate = NaT(sz);
            t.MergeAge = nan(sz);
            t.MergePtGender = cell(sz);
            t.MergePtEducat = nan(sz);
            t.MergePtEthCat = cell(sz);
            t.MergePtRacCat = cell(sz);
            t.MergeApoE4 = nan(sz);
            t.MergeFdg = nan(sz);
            t.MergePib = nan(sz);
            t.MergeAv45 = nan(sz);
            t.MergeAbeta = nan(sz);
            t.MergeTau = nan(sz);
            t.MergePTau = nan(sz);
            t.MergeMmse = nan(sz);
            t.MergeDx = cell(sz);
            t.MergeExamDateBl = NaT(sz);
            t.MergeFdgBl = nan(sz);
            t.MergePibBl = nan(sz);
            t.MergeAv45Bl = nan(sz);
            t.MergeAbetaBl = nan(sz);
            t.MergeTauBl = nan(sz);
            t.MergePTauBl = nan(sz);
            t.MergeMmseBl = nan(sz);
            t.MergeDxBl = cell(sz);
            t.MergeYearsBl = nan(sz);
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
                    merge_s1 = this.table_merge(strcmp(this.table_merge.PTID, s{1}), ':');
                    if isempty(merge_s1)
                        continue
                    end
                    rid = merge_s1.RID(1);
                    assert(all(rid == merge_s1.RID))
                    cdr_s1 = this.table_cdr(this.table_cdr.RID == rid, ':');
    
                    for acqdi = 1:size(t_s1, 1)
                        acqdate = t_s1.AcqDate(acqdi);

                        [~,merge_near] = min(abs(acqdate - merge_s1.EXAMDATE));
                        t_s1(acqdi, 'MergeVisCode') = merge_s1(merge_near, 'VISCODE');
                        t_s1(acqdi, 'MergeExamDate') = merge_s1(merge_near, 'EXAMDATE');
                        t_s1(acqdi, 'MergeAge') = merge_s1(merge_near, 'AGE');
                        t_s1(acqdi, 'MergePtGender') = merge_s1(merge_near, 'PTGENDER');
                        t_s1(acqdi, 'MergePtEducat') = merge_s1(merge_near, 'PTEDUCAT');
                        t_s1(acqdi, 'MergePtEthCat') = merge_s1(merge_near, 'PTETHCAT');
                        t_s1(acqdi, 'MergePtRacCat') = merge_s1(merge_near, 'PTRACCAT');
                        t_s1(acqdi, 'MergeApoE4') = merge_s1(merge_near, 'APOE4');
                        t_s1(acqdi, 'MergeFdg') = merge_s1(merge_near, 'FDG');
                        t_s1(acqdi, 'MergePib') = merge_s1(merge_near, 'PIB');
                        t_s1(acqdi, 'MergeAv45') = merge_s1(merge_near, 'AV45');
                        t_s1(acqdi, 'MergeAbeta') = merge_s1(merge_near, 'ABETA');
                        t_s1(acqdi, 'MergeTau') = merge_s1(merge_near, 'TAU');
                        t_s1(acqdi, 'MergePTau') = merge_s1(merge_near, 'PTAU');
                        t_s1(acqdi, 'MergeMmse') = merge_s1(merge_near, 'MMSE');
                        t_s1(acqdi, 'MergeDx') = merge_s1(merge_near, 'DX');
                        t_s1(acqdi, 'MergeExamDateBl') = merge_s1(merge_near, 'EXAMDATE_bl');
                        t_s1(acqdi, 'MergeFdgBl') = merge_s1(merge_near, 'FDG_bl');
                        t_s1(acqdi, 'MergePibBl') = merge_s1(merge_near, 'PIB_bl');
                        t_s1(acqdi, 'MergeAv45Bl') = merge_s1(merge_near, 'AV45_bl');
                        t_s1(acqdi, 'MergeAbetaBl') = merge_s1(merge_near, 'ABETA_bl');
                        t_s1(acqdi, 'MergeTauBl') = merge_s1(merge_near, 'TAU_bl');
                        t_s1(acqdi, 'MergePTauBl') = merge_s1(merge_near, 'PTAU_bl');
                        t_s1(acqdi, 'MergeMmseBl') = merge_s1(merge_near, 'MMSE_bl');
                        t_s1(acqdi, 'MergeDxBl') = merge_s1(merge_near, 'DX_bl');
                        t_s1(acqdi, 'MergeYearsBl') = merge_s1(merge_near, 'Years_bl');

                        [~,cdr_near] = min(abs(acqdate - cdr_s1.USERDATE));
                        t_s1(acqdi, 'VISCODE2') = cdr_s1{cdr_near, 'VISCODE2'};
                        t_s1{acqdi, 'USERDATE'} = cdr_s1{cdr_near, 'USERDATE'};
                        t_s1{acqdi, 'EXAMDATE'} = cdr_s1{cdr_near, 'EXAMDATE'};
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
                    end
                    
                    t(strcmp(t.Subject, s{1}), :) = t_s1;
                catch ME
                    handwarning(ME)
                end
            end
            %t = t(~isnan(t.RID), :);
        end
        function t = cdrRevisions(~, t)
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
