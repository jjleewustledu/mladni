classdef FDGDemographics < handle
    %% 
    %  
    %  Created 18-May-2022 20:19:23 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
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
                    j = jsondecode(fileread(jfile));
                    [~,id] = fileparts(j.ADNI_INFO.OriginalPath);
                    idids = [idids; id]; %#ok<AGROW> 
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
            icvs = mladni.FDGDemographics.find_icv(sarr);
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
                icv_ = cellfun(@(x) mladni.FDGDemographics.find_icv(x), obj, 'UniformOutput', false);
                icvs = cell2mat(icv_);
                icvs = ascol(icvs);
                return
            end
            assert(istext(obj));
            re = regexp(obj, '"image_mass":\s*(?<image_mass>\S+)', 'names');
            if isempty(re)
                icvs = 0;
                return
            end
            if length(re) > 1
                re = re(end);
            end
            icvs = str2double(re.image_mass);
            icvs = icvs*1.2; % ADNI has T1w voxels ~ [1.2 1 1]
        end
    end

    properties (Dependent)
        adni_demographics
    end
    
    methods

        %% GET
        
        function g = get.adni_demographics(this)
            g = this.adni_demo_;
        end
        
        %%
        
        function globbed_dx_csv = rawdata_pet_filename_dx(this, varargin)
            %% Builds and writes a subtable of rawdata PET filenames for subjects with a specified Merge Dx code.
            %  Args:
            %      globbed_csv (file):  default is /path/to/globbed.csv
            %      merge_dx (text):  'CN'|'Dementia'|'MCI' determined from AdniDemographics.
            
            globbed_csv = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.csv');
            
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
            
            globbed_csv = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.csv');
            
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
            
            globbed_csv = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.csv');
            
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
        function t = table_covariates(this, varargin)
            %% Saves and eturns table with ICV, PVE1 and all components in last column.
            %  Saves separate tables for each component.
            %
            %  Args:
            %      fn_csv (required file): csv containing single fqfns & component averages of all imaging files,
            %                              commonly 'component_weighted_average.csv'.
            %      trap_outliers (logical): removes scans with PVE1 > mean(PVE1) + 5*std(PVE1).
            %      save_1comp (logical): call this.table_covariates_1comp() for all components.

            ip = inputParser;
            addRequired(ip, 'fn_csv', @isfile);
            addParameter(ip, 'trap_outliers', true, @islogical);
            addParameter(ip, 'save_1comp', false, @islogical);
            parse(ip, varargin{:});
            ipr = ip.Results;

            % from ADNIDemographics
            t = table_fdg1(this);   
            t = removevars(t, ...
                {'Subject', 'RID', 'Sex', 'Age', 'Visit', 'Modality', 'Description', 'Type', 'AcqDate', 'Format', 'Downloaded', 'Group', ...
                 'MergeVisCode', 'MergeFdg', 'MergePib', 'MergeAv45', 'MergeAbeta', 'MergeTau', 'MergePTau', ...
                 'MergeExamDateBl', 'MergeFdgBl', 'MergePibBl', 'MergeAv45Bl', ...
                 'MergeAbetaBl', 'MergeTauBl', 'MergePTauBl', 'MergeMmseBl', 'MergeDxBl', 'MergeYearsBl', ...
                 'PonsVermis', ...
                 'VISCODE2', 'USERDATE', 'ID'});

            % ImageDataID
            [idids,comps,filelist] = mladni.FDGDemographics.csv_to_imagedataIDs(ipr.fn_csv);
            [idids,index] = sortrows(idids);
            Components = comps(index,:);
            Filelist = filelist(index);
            numComp = size(Components, 2);

            % ICV ~ intracranial mass
            icvs = mladni.FDGDemographics.csv_to_icvs(ipr.fn_csv);
            ICV = icvs(index);
            assert(length(Components) == length(ICV));

            % PVE1 ~ gray matter mass
            pve1s = mladni.FDGDemographics.csv_to_pve1(ipr.fn_csv);
            PVE1 = pve1s(index);
            assert(length(Components) == length(PVE1));

            % opportunistically trapping outliers
            if ipr.trap_outliers
                outliers = PVE1 < (mean(PVE1) - 5*std(PVE1));
                if ~isempty(Filelist(outliers))
                    for fn = Filelist(outliers)
                        fprintf('FDGDeomgraphics.table_covariates(): trapped outlier %s\n', fn{1});
                    end
                end
            else
                outliers = false(size(PVE1));
            end

            % assemble & save final table
            t = t(matches(t.ImageDataID, idids), :);
            t = sortrows(t, 1); % ImageDataID
            assert(all(string(t.ImageDataID) == string(idids)))
            t = addvars(t, Filelist, 'Before', 1);
            t = addvars(t, ICV, 'After', 'SITEID');
            t = addvars(t, PVE1, 'After', 'ICV');
            t = addvars(t, Components, 'After', 'PVE1');
            t = t(~outliers, :);
            t = t(~isnan(t.MergeAge),:);
            % t.ImageDataID = [];

            t = t(t.CDGLOBAL ~= -1, :);

            this.table_covariates_cache_ = t;
            save('mladni_FDGDemographics_table_covariates.mat', 't');
            writetable(t, 'mladni_FDGDemographics_table_covariates.csv');

            % save separate tables for each component
            if ipr.save_1comp
                for idx = 1:numComp
                    t1 = this.table_covariates_1comp(idx);
                    save(sprintf('mladni_FDGDemographics_table_covariates_1comp%i.mat', idx), 't1');
                end
            end
        end
        function t = table_covariates_1comp(this, idx)
            %% returns table with indexed component in last column.

            t = this.table_covariates_cache_;
            c = t.Components;
            c = c(:, idx);
            t.Components = c;
        end
        function t = table_fdg1(this)
            if ~isempty(this.table_fdg1_)
                t = this.table_fdg1_;
                return
            end
            this.table_fdg1_ = this.adni_demo_.table_fdg1();
            t = this.table_fdg1_;
        end
        
        function this = FDGDemographics(varargin)
            %% FDGDEMOGRAPHICS 
    
            ip = inputParser;
            addParameter(ip, 'table_covariates', []);
            parse(ip, varargin{:})
            ipr = ip.Results;

            this.adni_demo_ = mladni.AdniDemographics();
            this.table_covariates_cache_ = ipr.table_covariates;
        end
    end
    
    %% PRIVATE
    
    properties (Access = private)
        adni_demo_
        table_covariates_cache_
        table_fdg1_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
