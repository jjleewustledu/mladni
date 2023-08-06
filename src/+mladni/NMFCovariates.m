classdef NMFCovariates < handle
    %% Supports ADNI demographic data and other meta-data approximately at granularity of single scan objects.
    %  
    %  Created 18-May-2022 20:19:23 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Constant)
        EXCLUSIONS = "" % "sub-022S0096_ses-20060228101016_trc-FDG_proc-CASU-ponsvermis_orient-rpi_pet_on_T1w_Warped_dlicv.nii.gz"
    end

    properties (Dependent)
        componentDir
        covariates_file
        covariates_1stscan_file
        inFiles
        selectedNumBases
        study_design
        targetDatasetDir
        weightedAverFilename
        X2mat
    end
    
    methods % GET
        function g = get.componentDir(this)
            g = fullfile(this.targetDatasetDir, sprintf('NumBases%i', this.selectedNumBases), 'components');
            ensuredir(g);
        end
        function g = get.covariates_file(this)
            g = fullfile(this.componentDir, sprintf("NMFCovariates_table_covariates_%s.mat", this.study_design));
        end
        function g = get.covariates_1stscan_file(this)
            g = fullfile(this.componentDir, sprintf("NMFCovariates_table_covariates_1stscan_%s.mat", this.study_design));
        end
        function g = get.inFiles(this)
            g = fullfile(this.componentDir, sprintf("%s_%s.csv", stackstr(), this.study_design));
            if ~isfile(g)
                t = table(this.table_fdg4.Filelist);
                writetable(t, g, WriteVariableNames=false)
            end
        end
        function g = get.selectedNumBases(this)
            g = this.selectedNumBases_;
        end
        function g = get.study_design(this)
            g = this.study_design_;
        end
        function g = get.targetDatasetDir(~)
            g = fullfile(getenv("ADNI_HOME"), "NMF_FDG", "baseline_cn");
        end
        function g = get.weightedAverFilename(this)
            g = fullfile(this.componentDir, ...
                sprintf("component_weighted_average_%s.csv", this.study_design));
        end
        function g = get.X2mat(this)
            g = fullfile(this.componentDir, sprintf("X2_%s.mat", this.study_design));
        end
    end

    methods
        function this = NMFCovariates(opts)
            %% NMFCovariates 
            %  Args:
            %    opts.table_covariates = []: to provision from memory
            %    opts.study_design = "cross-sectional"            
    
            arguments
                opts.table_covariates = []
                opts.selectedNumBases {mustBeScalarOrEmpty} = 22
                opts.study_design = "longitudinal"
            end
            this.selectedNumBases_ = opts.selectedNumBases;
            this.study_design_ = opts.study_design;
            this.adni_demo_ = mladni.AdniDemographics(study_design=opts.study_design);
            this.table_covariates_cache_ = opts.table_covariates;
            %assert(this.components_are_available)
        end   
        function globbed_dx_csv = rawdata_pet_filename_dx(this, varargin)
            %% Builds and writes a subtable of rawdata PET filenames for subjects with a specified Merge Dx code.
            %  Args:
            %      globbed_csv (file):  default is /path/to/globbed.csv
            %      merge_dx (text):  'CN'|'Dementia'|'MCI' determined from AdniDemographics.
            
            globbed_csv = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'mladni_FDG_batch_globbed.csv');
            
            ip = inputParser;
            addParameter(ip, 'globbed_csv', globbed_csv, @isfile);
            addParameter(ip, 'merge_dx', 'CN', @(x) any(contains({'CN', 'Dementia', 'MCI'}, x)))
            addParameter(ip, 'description', 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;           
            globbed_tbl = readtable(ipr.globbed_csv, 'ReadVariableNames', false, 'Delimiter', ' '); % 3735 x 1 table
            
            globbed_dx_csv = strcat(myfileprefix(ipr.globbed_csv), '_', lower(ipr.merge_dx), '.csv');
            
            %% select dx from table_fdg1
            
            fdg1 = this.adni_demo_.table_fdg1(); % 14358 x 53 table
            dx = fdg1.MergeDx; % 15358 x 1 cell
            desc = fdg1.Description; % 15358 x 1 cell
            select1 = cell2mat( ...
                cellfun(@(x) ischar(x) && strcmp(x, ipr.description), desc, ...
                'UniformOutput', false));
            select2 = cell2mat( ...
                cellfun(@(x) ischar(x) && strcmp(x, ipr.merge_dx), dx, ...
                'UniformOutput', false));
            select = select1 & select2;
            sub = fdg1.Subject(select);
            sub = cellfun(@(x) strrep(x, '_', ''), sub, 'UniformOutput', false);
            ses = fdg1.AcqDate(select);
            ses = cellstr(datestr(ses, 'yyyymmdd'));
            
            %% select sub and ses from globbed_tbl
            
            selected_files = {};
            files = globbed_tbl.Var1;
            for rowi = 1:length(files)
                re = regexp(files{rowi}, '\S+/sub-(?<sub>\d{3}S\d{4})/ses-(?<ses>\d{8})/\S+', 'names');
                if sum( contains(sub, re.sub) & contains(ses, re.ses) ) > 0
                    selected_files = [selected_files; files{rowi}]; %#ok<AGROW>
                end
            end
            
            %% write requested subtable
            
            tbl = table(selected_files);
            tbl.Properties.VariableNames = ...
                cellfun(@(x) strcat(x, '_', lower(ipr.merge_dx)), globbed_tbl.Properties.VariableNames, ...
                'UniformOutput', false);
            writetable(tbl, globbed_dx_csv, 'WriteVariableNames', false);
        end        
        function [globbed_amypos_csv,globbed_amyneg_csv] = rawdata_pet_filename_amyloid(this, varargin)
            %% Builds and writes subtables of rawdata PET filenames for subjects with and without amyloid;
            %  positivity threshed per 
            %  https://adni.bitbucket.io/reference/docs/UCBERKELEYAV45/UCBERKELEY_AV45_Methods_04.25.2022.pdf
            %  Args:
            %      globbed_csv (file):  default is /path/to/globbed.csv
            %      description (text):  e.g., 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution'.
            %  Returns:
            %      globbed_amypos_csv:  filename.
            %      globbed_amyneg_csv:  filename.
            
            globbed_csv = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'mladni_FDG_batch_globbed.csv');
            
            ip = inputParser;
            addParameter(ip, 'globbed_csv', globbed_csv, @isfile);
            addParameter(ip, 'description', 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;           
            globbed_tbl = readtable(ipr.globbed_csv, 'ReadVariableNames', false, 'Delimiter', ' '); % 3735 x 1 table
            
            globbed_amypos_csv = strcat(myfileprefix(ipr.globbed_csv), '_amypos.csv');
            globbed_amyneg_csv = strcat(myfileprefix(ipr.globbed_csv), '_amyneg.csv');
            
            %% select amyloid status using ADNI/studydata/ucberkeleyav45_skinny.csv
            
            fdg1 = this.adni_demo_.table_fdg1(); % 14358 x 53 table
            sub = cellfun(@(x) strrep(x, '_', ''), fdg1.Subject, 'UniformOutput', false);
            ses = cellstr(datestr(fdg1.AcqDate, 'yyyymmdd'));

            desc = fdg1.Description; % 15358 x 1 cell
            select_desc = cell2mat( ...
                cellfun(@(x) ischar(x) && strcmp(x, ipr.description), desc, ...
                'UniformOutput', false));

            amyloid = fdg1.AmyloidStatus; % 15358 x 1 logical
            amyloid(isnan(amyloid)) = -1;
            select_amypos = (amyloid == 1) & select_desc;
            select_amyneg = (amyloid == 0) & select_desc;
            sub_amypos = sub(select_amypos);
            sub_amyneg = sub(select_amyneg);
            ses_amypos = ses(select_amypos);
            ses_amyneg = ses(select_amyneg);
            
            %% select sub and ses from globbed_tbl
            
            files_amypos = {};
            files_amyneg = {};
            files = globbed_tbl.Var1;
            for rowi = 1:length(files)
                re = regexp(files{rowi}, '\S+/sub-(?<sub>\d{3}S\d{4})/ses-(?<ses>\d{8})/\S+', 'names');
                if sum(contains(sub_amypos, re.sub) & contains(ses_amypos, re.ses)) > 0
                    files_amypos = [files_amypos; files{rowi}]; %#ok<AGROW>
                end
                if sum(contains(sub_amyneg, re.sub) & contains(ses_amyneg, re.ses)) > 0
                    files_amyneg = [files_amyneg; files{rowi}]; %#ok<AGROW>
                end
            end
            
            %% write requested subtables
            
            tbl_amypos = table(files_amypos);
            tbl_amyneg = table(files_amyneg);
            writetable(tbl_amypos, globbed_amypos_csv, 'WriteVariableNames', false);
            writetable(tbl_amyneg, globbed_amyneg_csv, 'WriteVariableNames', false);
        end
        function globbed_dx_csv = pve_filename_dx(this, varargin)
            %% Builds and writes a subtable of pve filenames for subjects with a specified Merge Dx code.
            %  Args:
            %      globbed_csv (file):  default is /path/to/globbed.csv
            %      merge_dx (text):  'CN'|'Dementia'|'MCI' determined from AdniDemographics.
            
            globbed_csv = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'mladni_FDG_batch_globbed.csv');
            
            ip = inputParser;
            addParameter(ip, 'globbed_csv', globbed_csv, @isfile);
            addParameter(ip, 'merge_dx', 'CN', @(x) any(contains({'CN', 'Dementia', 'MCI'}, x)))
            addParameter(ip, 'description', 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution', @istext)
            addParameter(ip, 'pve', 1, @isscalar)
            parse(ip, varargin{:});
            ipr = ip.Results;           
            globbed_tbl = readtable(ipr.globbed_csv, 'ReadVariableNames', false, 'Delimiter', ' '); % 3735 x 1 table
            
            globbed_dx_csv = sprintf('%s_pve%i_%s.csv', ...
                myfileprefix(ipr.globbed_csv), ipr.pve, lower(ipr.merge_dx));
            
            %% select dx from table_fdg1
            
            fdg1 = this.adni_demo_.table_fdg1(); % 14358 x 53 table
            dx = fdg1.MergeDx; % 15358 x 1 cell
            desc = fdg1.Description; % 15358 x 1 cell
            select1 = cell2mat( ...
                cellfun(@(x) ischar(x) && strcmp(x, ipr.description), desc, ...
                'UniformOutput', false));
            select2 = cell2mat( ...
                cellfun(@(x) ischar(x) && strcmp(x, ipr.merge_dx), dx, ...
                'UniformOutput', false));
            select = select1 & select2;
            sub = fdg1.Subject(select);
            sub = cellfun(@(x) strrep(x, '_', ''), sub, 'UniformOutput', false);
            ses = fdg1.AcqDate(select);
            ses = cellstr(datestr(ses, 'yyyymmdd'));
            
            %% select sub and ses from globbed_tbl
            
            selected_files = {};
            files = globbed_tbl.rawdata_pet_filename;
            for rowi = 1:length(files)
                re = regexp(files{rowi}, '\S+/sub-(?<sub>\d{3}S\d{4})/ses-(?<ses>\d{8})/\S+', 'names');
                if sum( contains(sub, re.sub) & contains(ses, re.ses) ) > 0
                    pth = strrep(myfileparts(files{rowi}), 'rawdata', 'derivatives');
                    g = glob(fullfile(pth, ...
                        sprintf('sub-%s_ses-*_*_orient-rpi_T1w_brain_pve_%i_detJ_Warped.nii.gz', re.sub, ipr.pve)));
                    if ~isempty(g)
                        selected_files = [selected_files; g{1}]; %#ok<AGROW>
                    end
                end
            end
            
            %% write requested subtable
            
            tbl = table(selected_files);
            tbl.Properties.VariableNames = {sprintf('filename_pve%i_%s', ipr.pve, lower(ipr.merge_dx))};
            writetable(tbl, globbed_dx_csv);
        end        
        
        %% tables

        function t = AddAcqDuration(~, t)
            %% AcqDuration ~ floating-point years since first scan
            %  Age ~ adjusted to precision of AcqDuration after first scan

            nats = nan(size(t.AcqDate));
            t = addvars(t, nats, NewVariableNames={'AcqDuration'}, After='AcqDate');
            subs = unique(t.Subject);
            for s = asrow(subs)
                sub_select = strcmp(t.Subject, s{1});
                u = t(sub_select, :);
                assert(all(u.AcqDate >= u.AcqDate(1)))
                acqdur = years(u.AcqDate - u.AcqDate(1));
                adjage = u.Age(1) + acqdur;
                t.AcqDuration(sub_select) = acqdur;
                t.Age(sub_select) = adjage;
            end
        end
        function t = addvars_by_filelist(this, t, t_, new_var_, opts)
            arguments
                this mladni.NMFCovariates
                t table % containing longer Filelist
                t_ table % containing shorter Filelist
                new_var_ double % matching t_
                opts.NewVariableNames cell 
            end
            assert(any(contains(t.Properties.VariableNames, "Filelist")), stackstr())
            assert(any(contains(t_.Properties.VariableNames, "Filelist")), stackstr())
            assert(size(t,1) >= size(t,2))

            select = contains(t.Filelist, t_.Filelist);
            new_var = nan(size(t,1), size(new_var_,2));
            new_var(select,:) = new_var_;
            t = addvars(t, new_var, NewVariableNames=opts.NewVariableNames);
        end
        function t = apply_table_qc(this, t)
            vns = t.Properties.VariableNames;
            
            % start with 3478 rows

            if any(contains(vns, 'CDGLOBAL'))
                t = t(t.CDGLOBAL ~= -1, :);
            end % 
            if any(contains(vns, 'Cohort'))
                t = t(t.Cohort ~= 'unknown', :);
            end % 

            if any(contains(vns, 'Components_1'))
                t = t(~any(isnan(t.Components_1),2), :);
            end % 
            if any(contains(vns, 'Age'))
                t = t(~isnan(t.Age), :);
            end % 
            if any(contains(vns, 'ApoE4'))
                t = t(~isnan(t.ApoE4), :);
            end % 
            
%             if any(contains(vns, 'Dlicv'))
%                 t = t(t.Dlicv > 0.5e6, :);
%             end % 3464 rows
%             if any(contains(vns, 'PVE1'))
%                 t = t(t.PVE1 > 0.5e6/3, :);
%             end % 3464 rows
%             if any(contains(vns, 'RegErr'))
%                 t = t(t.RegErr < 5, :); % see also mladni.FDG()
%             end % 3435 rows

            if ~isemptytext(this.EXCLUSIONS)
                for fidx = 1:length(this.EXCLUSIONS)
                    select = ~contains(t.Filelist, this.EXCLUSIONS(fidx));
                    t = t(select, :);
                end
            end
        end
        function t = table_covariates(this, opts)
            %% Saves and returns table useful for inferences using RStudio.
            %  Args:
            %      save_1comp (logical): call this.table_covariates_1comp() for all components.

            arguments
                this mladni.NMFCovariates
                opts.save_1comp logical = false                
            end

            if ~isempty(this.table_covariates_cache_)
                fprintf("%s: using cached in memory\n", stackstr())
                t = this.table_covariates_cache_;
                return
            end
            cache_file = this.covariates_file;
            if isfile(cache_file)
                fprintf("%s: using cached from filesystem\n", stackstr())
                ld = load(cache_file);
                this.table_covariates_cache_ = ld.t;
                t = this.table_covariates_cache_;
                return
            end

            % from ADNIDemographics ~ 3478 rows
            t = table_fdg4(this); 

            % Cohort ~ categorical
            Cohort = repmat({'unknown'}, size(t,1), 1);
            Cohort(t.CDGLOBAL == 0 & t.AmyloidStatusLong == 0) = {'CDR=0,amy-'};
            Cohort(t.CDGLOBAL == 0 & t.AmyloidStatusLong == 1) = {'CDR=0,amy+'};
            Cohort(t.CDGLOBAL == 0.5 & t.AmyloidStatusLong == 1) = {'CDR=0.5,amy+'};
            Cohort(t.CDGLOBAL > 0.5 & t.AmyloidStatusLong == 1) = {'CDR>0.5,amy+'};
            Cohort(t.CDGLOBAL > 0 & t.AmyloidStatusLong == 0) = {'CDR>0,amy-'};
            Cohort = categorical(Cohort);
            t = addvars(t, Cohort, NewVariableNames={'Cohort'});

            % DLICV ~ 3478 rows, 12 nans
            t_ = this.table_dlicv;
            t = addvars(t, t_.Dlicv, NewVariableNames={'Dlicv'});

            % PVE1 ~ 3478 rows, 8 nans
            t_ = this.table_pve1;
            t = addvars(t, t_.PVE1, NewVariableNames={'PVE1'});

            % RegErr ~ 3478 rows, 38 nans
            t_ = this.table_regerr;
            t = addvars(t, t_.RegErr, NewVariableNames={'RegErr'});

            % Components ~ 3470 rows -> 3478 rows, 8 nans
            % implicitly marks empty Filelist for exclusion by apply_table_qc()
            t_ = this.table_selectedComponentWeightedAverageNIFTI;
            t = this.addvars_by_filelist(t, t_, t_.Components, NewVariableNames={'Components'});
            t = splitvars(t, 'Components');

            % apply table qc
            t = this.apply_table_qc(t);

            % sort rows by Subject, then AcqDate
            t = sortrows(t, ["Subject", "AcqDate"]);

            % AcqDuration ~ floating-point years since first scan
            t = this.AddAcqDuration(t);

            % store cache, & save/write table
            this.table_covariates_cache_ = t;
            save(cache_file, 't');
            writetable(t, strrep(cache_file, ".mat", ".csv"));

            % save separate tables for each component
            if opts.save_1comp
                for idx = 1:this.selectedNumBases
                    t1 = this.table_covariates_1comp(idx);
                    save(strrep(cache_file, ".mat", "_1comp.mat"), 't1');
                    writetable(t1, strrep(cache_file, ".mat", "_1comp.csv"));
                end
            end
        end
        function t = table_covariates_1stscan(this)
            t = this.adni_demo_.table_firstscan(this.table_covariates());

            cache_file = this.covariates_1stscan_file;
            save(cache_file, 't');
            writetable(t, strrep(cache_file, ".mat", ".csv"));
        end
        function t = table_covariates_1comp(this, idx)
            %% returns table with indexed component in last column.

            t = this.table_covariates_cache_;
            c = t.Components;
            c = c(:, idx);
            t.Components = c;
        end
        function t = table_fdg4(this)
            if ~isempty(this.table_fdg4_)
                t = this.table_fdg4_;
                return
            end
            this.table_fdg4_ = this.adni_demo_.table_fdg4();
            t = this.table_fdg4_;
        end
        function t = table_filelist(this)
            t = table(this.table_fdg4.Filelist);
            t.Properties.VariableNames = {'Filelist'};
        end
        function t = table_dlicv(this)
            %% Finds dlicv json string; estimate icv using this.find_icv().
            %  Returns:
            %      t: table with variables Filelist and DLICV.

            Filelist = this.table_fdg4.Filelist;
            Dlicv = nan(size(Filelist));
            for fidx = 1:length(Filelist)
                item = Filelist{fidx};
                if isempty(item)
                    continue; 
                end
                try
                    if contains(item, ',')
                        ss = strsplit(item, ',');
                        item = ss{1};
                    end
                    if ~contains(item, 'orient-rpi_T1w_dlicv.nii.gz')
                        pth = fileparts(item);
                        g = glob(fullfile(pth, '*orient-rpi_T1w_dlicv.nii.gz'));
                        assert(~isempty(g))
                        item = g{1};
                    end
                    jfile = strrep(item, '.nii.gz', '.json');
                    Dlicv(fidx) = mladni.NMFCovariates.find_icv(fileread(jfile));
                catch ME
                    handwarning(ME);
                end
            end
            t = table(Filelist, Dlicv);
        end
        function t = table_imagedataIDs(this)
            %% Finds imagedataIDs;
            %  Returns:
            %      t: table with variables Filelist and ImageDataID.

            Filelist = this.table_fdg4.Filelist;
            ImageDataID = cell(size(Filelist));
            for fidx = 1:length(Filelist)
                ImageDataID{fidx} = 'unknown';
                item = Filelist{fidx};
                if isempty(item)
                    continue; 
                end
                try
                    itemj = '';
                    if contains(item, ',')
                        ss = strsplit(item, ',');
                        item = ss{1};
                        itemj = item;
                    end
                    if ~contains(item, 'orient-rpi_pet')
                        pth = fileparts(item);
                        g = glob(fullfile(pth, '*orient-rpi_pet.json'));
                        assert(~isempty(g))
                        itemj = g{1};
                    end
                    if isempty(itemj)
                        itemj = item;
                    end
                    jfile = strrep(itemj, '.nii.gz', '.json');
                    if ~isfile(jfile)
                        % KLUDGE
                        itemj = strrep(itemj, '_dlicv', '');
                        jfile = strrep(itemj, '.nii.gz', '.json');
                    end
                    if ~isfile(jfile)
                        % KLUDGE
                        itemj = strrep(itemj, '_Warped', '');
                        jfile = strrep(itemj, '.nii.gz', '.json');
                    end
                    if ~isfile(jfile)
                        % KLUDGE
                        itemj = strrep(itemj, '-ponsvermis', '');
                        jfile = strrep(itemj, '.nii.gz', '.json');
                    end
                    j = jsondecode(fileread(jfile));
                    if endsWith(j.ADNI_INFO.OriginalPath, filesep)
                        opath = convertStringsToChars(j.ADNI_INFO.OriginalPath);
                        j.ADNI_INFO.OriginalPath = opath(1:end-1);
                    end
                    [~,id] = fileparts(j.ADNI_INFO.OriginalPath);
                    ImageDataID{fidx} = id;
                catch ME
                    handwarning(ME);
                end
            end
            t = table(Filelist, ImageDataID);
        end
        function t = table_pve1(this)
            %% Finds image_mass json string; estimate mass of pve1.
            %  Returns:
            %      t: table with variables Filelist and PVE1, grey matter.

            Filelist = this.table_fdg4.Filelist;
            PVE1 = nan(size(Filelist));
            for fidx = 1:length(Filelist)
                item = Filelist{fidx};
                if isempty(item)
                    continue; 
                end
                try
                    if contains(item, ',')
                        ss = strsplit(item, ',');
                        item = ss{1};
                    end
                    if ~contains(item, 'orient-rpi_T1w_brain_pve_1.nii.gz')
                        pth = fileparts(item);
                        g = glob(fullfile(pth, '*orient-rpi_T1w_brain_pve_1.nii.gz'));
                        assert(~isempty(g))
                        item = g{1};
                    end
                    jfile = strrep(item, '.nii.gz', '.json');
                    s = readlines(jfile);
                    j = jsondecode([s{:}]);
                    PVE1(fidx) = j.FDG_fast.image_mass*1.2; % ADNI MPRAGE has voxel mmppix ~ [1.2 1 1]
                catch ME
                    handwarning(ME);
                end
            end
            t = table(Filelist, PVE1);
        end
        function t = table_regerr(this)
            Filelist = this.table_fdg4.Filelist;
            RegErr = nan(size(Filelist));
            for fidx = 1:length(Filelist)
                item = Filelist{fidx};
                if isempty(item)
                    continue; 
                end
                try
                    if contains(item, ',')
                        ss = strsplit(item, ',');
                        item = ss{1};
                    end
                    if ~contains(item, 'orient-rpi_pet_on_T1w.nii.gz')
                        pth = fileparts(item);
                        g = glob(fullfile(pth, '*orient-rpi_pet_on_T1w.nii.gz'));
                        assert(~isempty(g))
                        item = g{1};
                    end
                    jfile = strrep(item, '.nii.gz', '.json');
                    RegErr(fidx) = mladni.NMFCovariates.find_regerr(fileread(jfile));
                catch ME
                    handwarning(ME);
                end
            end
            t = table(Filelist, RegErr);
        end
        function t = table_selectedComponentWeightedAverageNIFTI(this)
            %% See also mladni.NMF.calculateSelectedComponentWeightedAverageNIFTI(),
            %  which writes component_weighted_average_<study_design>.csv.
            %  t.Filelist will exclude registration failures.

            arguments
                this mladni.NMFCovariates
            end

            % prep param struct
            param.isList = 1 ;
            param.downSample = 1;
            param.smooth = 0;
            param.mask = fullfile(getenv("ADNI_HOME"), "VolBin", "mask.nii.gz");
            param.permute = 0;
            param.numBase = this.selectedNumBases;
            param.componentDir = this.componentDir;
            assert(isfile(param.mask))

            % ensure X2.mat
            if isfile(this.X2mat)
                load(this.X2mat);
            else
                data = mladni.NMF.loadData(this.inFiles, param, []);
                meanX = mean(data.X, 2);
                X = data.X(meanX > 0, :); 
                clear data;
                save(this.X2mat, "X", "meanX");
            end

            % generate table t; write table
            if ~isfile(this.weightedAverFilename)
                mladni.NMF.calculateSelectedComponentWeightedAverageNIFTI( ...
                    this.inFiles, this.targetDatasetDir, this.selectedNumBases, this.weightedAverFilename);
            end
            t = readtable(this.weightedAverFilename);
            comp_vars = cellfun(@(x) sprintf('P%i', x), num2cell(1:this.selectedNumBases), UniformOutput=false);
            t.Properties.VariableNames = ['Filelist', comp_vars];
            t = mergevars(t, comp_vars, NewVariableName="Components");
        end

        %% tables of diagnostic subgroups

        function t = table_all(this, cs, T_name)
            %% N -> 3415 longitudinal, 1529 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = strcmp(this.study_design, "cross-sectional")
                T_name {mustBeTextScalar} = 'table_covariates'
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__;
            else
                t = this.table_firstscan(fdg__);
            end
        end
        function t = table_cn(this, cs, T_name, bl_1st)
            %% N -> 474 longitudinal, 264 cross-sectional.
            
            arguments
                this mladni.NMFCovariates
                cs logical = strcmp(this.study_design, "cross-sectional")
                T_name {mustBeTextScalar} = 'table_covariates'
                bl_1st logical = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusLong == 0, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusLong == 0, :);
            else
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusLong == 0, :);
                t = this.table_firstscan(t);
            end
        end
        function t = table_preclinical(this, cs, T_name, bl_1st)
            %% N -> 174 longitudinal, 137 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = strcmp(this.study_design, "cross-sectional")
                T_name {mustBeTextScalar} = 'table_covariates'
                bl_1st logical = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusLong == 1, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusLong == 1, :);
            else
                t = fdg__(fdg__.CDGLOBAL == 0 & fdg__.AmyloidStatusLong == 1, :);
                t = this.table_firstscan(t);
            end
        end        
        function t = table_cdr_0p5_apos(this, cs, T_name, bl_1st)
            %% N -> 564 longitudinal, 455 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = strcmp(this.study_design, "cross-sectional")
                T_name {mustBeTextScalar} = 'table_covariates'
                bl_1st logical = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL == 0.5 & fdg__.AmyloidStatusLong == 1, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL == 0.5 & fdg__.AmyloidStatusLong == 1, :);
            else
                t = fdg__(fdg__.CDGLOBAL == 0.5 & fdg__.AmyloidStatusLong == 1, :);
                t = this.table_firstscan(t);
            end
        end
        function t = table_cdr_gt_0p5_apos(this, cs, T_name, bl_1st)
            %% N -> 174 longitudinal, 150 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = strcmp(this.study_design, "cross-sectional")
                T_name {mustBeTextScalar} = 'table_covariates'
                bl_1st logical = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL > 0.5 & fdg__.AmyloidStatusLong == 1, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL > 0.5 & fdg__.AmyloidStatusLong == 1, :);
            else
                t = fdg__(fdg__.CDGLOBAL > 0.5 & fdg__.AmyloidStatusLong == 1, :);
                t = this.table_firstscan(t);
            end            
        end
        function t = table_cdr_gt_0_aneg(this, cs, T_name, bl_1st)
            %% N -> 557 longitudinal, 334 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = strcmp(this.study_design, "cross-sectional")
                T_name {mustBeTextScalar} = 'table_covariates'
                bl_1st logical = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL >= 0.5 & fdg__.AmyloidStatusLong == 0, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL >= 0.5 & fdg__.AmyloidStatusLong == 0, :);
            else
                t = fdg__(fdg__.CDGLOBAL >= 0.5 & fdg__.AmyloidStatusLong == 0, :);
                t = this.table_firstscan(t);
            end            
        end 
        function t = table_firstscan(this, varargin)
            t = this.adni_demo_.table_firstscan(varargin{:});
        end
    end

    methods (Static)
        function create_tables_for_R()
        end
        function tf = components_are_available()
            tf = mladni.NMF.components_are_available();
        end
        function [idids,comps,filelist] = csv_to_imagedataIDs(fn_csv)
            %% finds imagedataIDs from csv;
            %  also finds numerical components if fn_csv is, e.g., "component_weighted_average.csv"
            %  Params:
            %      fn_csv (file): containing fqfn for NIfTI
            %
            %  Returns:
            %      idids: cell array of char
            %      comps: numeric component weighted averages
            %      filelist: cell array of filenames from from fn_csv

            lines = readlines(fn_csv);
            idids = {};
            filelist = {};
            comps = [];
            for f = asrow(lines)
                if isempty(f{1}); continue; end
                try
                    item = f{1};
                    itemj = '';
                    if contains(item, ',')
                        ss = strsplit(item, ',');
                        item = ss{1};
                        itemj = item;
                        comps = [comps; str2double(ss(2:end))]; %#ok<AGROW> 
                    end
                    if ~contains(item, 'orient-rpi_pet')
                        pth = fileparts(item);
                        g = glob(fullfile(pth, '*orient-rpi_pet.json'));
                        assert(~isempty(g))
                        itemj = g{1};
                    end
                    if isempty(itemj)
                        itemj = item;
                    end
                    jfile = strrep(itemj, '.nii.gz', '.json');
                    if ~isfile(jfile)
                        % KLUDGE
                        itemj = strrep(itemj, '_dlicv', '');
                        jfile = strrep(itemj, '.nii.gz', '.json');
                    end
                    if ~isfile(jfile)
                        % KLUDGE
                        itemj = strrep(itemj, '_Warped', '');
                        jfile = strrep(itemj, '.nii.gz', '.json');
                    end
                    if ~isfile(jfile)
                        % KLUDGE
                        itemj = strrep(itemj, '-ponsvermis', '');
                        jfile = strrep(itemj, '.nii.gz', '.json');
                    end
                    try
                        j = jsondecode(fileread(jfile));
                        if endsWith(j.ADNI_INFO.OriginalPath, filesep)
                            opath = convertStringsToChars(j.ADNI_INFO.OriginalPath);
                            j.ADNI_INFO.OriginalPath = opath(1:end-1);
                        end
                        [~,id] = fileparts(j.ADNI_INFO.OriginalPath);
                        idids = [idids; id]; %#ok<AGROW> 
                    catch
                        idids = [idids; 'unknown']; %#ok<AGROW> 
                    end
                    filelist = [filelist; item]; %#ok<AGROW> 
                catch ME
                    handwarning(ME);
                    if ~isempty(comps)
                        comps = comps(1:end-1,:);
                    end
                end
            end
        end
        function icvs = csv_to_icvs(fn_csv)
            %% finds dlicv json string; estimate icv from image_mass
            %  Params:
            %      fn_csv (file): containing fqfn for NIfTI
            %
            %  Returns:
            %      icvs: vector of icvs

            idx = 1;
            filelist = readlines(fn_csv);
            for f = asrow(filelist)
                if isempty(f{1}); continue; end
                try
                    item = f{1};
                    if contains(item, ',')
                        ss = strsplit(item, ',');
                        item = ss{1};
                    end
                    if ~contains(item, 'orient-rpi_T1w_dlicv.nii.gz')
                        pth = fileparts(item);
                        g = glob(fullfile(pth, '*orient-rpi_T1w_dlicv.nii.gz'));
                        assert(~isempty(g))
                        item = g{1};
                    end
                    jfile = strrep(item, '.nii.gz', '.json');
                    sarr{idx} = fileread(jfile); %#ok<AGROW> 
                catch ME
                    handwarning(ME);
                    sarr{idx} = ""; %#ok<AGROW> 
                end
                idx = idx + 1;
            end
            icvs = mladni.NMFCovariates.find_icv(sarr);
        end
        function pve1s = csv_to_pve1(fn_csv)
            %% finds mass of pve1 from csv
            %  Params:
            %      fn_csv (file): containing fqfn for NIfTI
            %
            %  Returns:
            %      pve1s: vector of mass of pve1

            idx = 1;
            filelist = readlines(fn_csv);
            for f = asrow(filelist)
                if isempty(f{1}); continue; end
                try
                    item = f{1};
                    if contains(item, ',')
                        ss = strsplit(item, ',');
                        item = ss{1};
                    end
                    if ~contains(item, 'orient-rpi_T1w_brain_pve_1.nii.gz')
                        pth = fileparts(item);
                        g = glob(fullfile(pth, '*orient-rpi_T1w_brain_pve_1.nii.gz'));
                        assert(~isempty(g))
                        item = g{1};
                    end
                    jfile = strrep(item, '.nii.gz', '.json');
                    s = readlines(jfile);
                    j = jsondecode([s{:}]);
                    pve1s(idx) = j.FDG_fast.image_mass*1.2; %#ok<AGROW> 
                catch ME
                    handwarning(ME);
                    pve1s(idx) = 0; %#ok<AGROW> 
                end
                idx = idx + 1;
            end
            pve1s = ascol(pve1s);
        end
        function icvs = find_icv(obj)
            if iscell(obj)
                icv_ = cellfun(@(x) mladni.NMFCovariates.find_icv(x), obj, 'UniformOutput', false);
                icvs = cell2mat(icv_);
                icvs = ascol(icvs);
                return
            end
            assert(istext(obj));
            re = regexp(obj, '"image_mass":\s*(?<image_mass>\S+)', 'names');
            if isempty(re)
                icvs = nan;
                return
            end
            if length(re) > 1
                re = re(end);
            end
            icvs = str2double(re.image_mass);
            icvs = icvs*1.2; % ADNI has T1w voxels ~ [1.2 1 1]
        end
        function errs = find_regerr(obj)
            if iscell(obj)
                err_ = cellfun(@(x) mladni.NMFCovariates.find_regerr(x), obj, 'UniformOutput', false);
                errs = cell2mat(err_);
                errs = ascol(errs);
                return
            end
            assert(istext(obj));
            re = regexp(obj, '"err":\s*"(?<err>\S+)"', 'names');
            if isempty(re)
                errs = nan;
                return
            end
            if length(re) > 1
                re = re(end);
            end
            errs = str2double(re.err);
        end
    end
    
    %% PRIVATE
    
    properties (Access = private)
        adni_demo_
        selectedNumBases_
        study_design_
        table_cohorts_
        table_covariates_cache_
        table_fdg4_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
