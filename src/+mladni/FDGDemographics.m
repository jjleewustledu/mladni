classdef FDGDemographics < handle
    %% 
    %  
    %  Created 18-May-2022 20:19:23 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
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
            globbed_tbl = readtable(ipr.globbed_csv, 'ReadVariableNames', true, 'Delimiter', ' '); % 3735 x 1 table
            
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
            files = globbed_tbl.rawdata_pet_filename;
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
            writetable(tbl, globbed_dx_csv);
        end        
        function globbed_amyloid_csv = rawdata_pet_filename_amyloid(this, varargin)
            %% Builds and writes a subtable of rawdata PET filenames for subjects with amyloid;
            %  positivity threshed per 
            %  https://adni.bitbucket.io/reference/docs/UCBERKELEYAV45/UCBERKELEY_AV45_Methods_04.25.2022.pdf
            %  Args:
            %      globbed_csv (file):  default is /path/to/globbed.csv
            
            globbed_csv = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.csv');
            
            ip = inputParser;
            addParameter(ip, 'globbed_csv', globbed_csv, @isfile);
            addParameter(ip, 'description', 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;           
            globbed_tbl = readtable(ipr.globbed_csv, 'ReadVariableNames', true, 'Delimiter', ' '); % 3735 x 1 table
            
            globbed_amyloid_csv = strcat(myfileprefix(ipr.globbed_csv), '_amyloid.csv');
            
            %% select amyloid status using ADNI/studydata/ucberkeleyav45_skinny.csv
            
            fdg1 = this.adni_demo_.table_fdg1(); % 14358 x 53 table
            amyloid = fdg1.AmyloidBinary; % 15358 x 1 logical
            desc = fdg1.Description; % 15358 x 1 cell
            select1 = cell2mat( ...
                cellfun(@(x) ischar(x) && strcmp(x, ipr.description), desc, ...
                'UniformOutput', false));
            select2 = cell2mat( ...
                cellfun(@(x) islogical(x) && x, amyloid, ...
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
                    selected_files = [selected_files; files{rowi}]; %#ok<AGROW>
                end
            end
            
            %% write requested subtable
            
            tbl = table(selected_files);
            tbl.Properties.VariableNames = ...
                cellfun(@(x) strcat(x, '_amyloid', globbed_tbl.Properties.VariableNames), ...
                'UniformOutput', false);
            writetable(tbl, globbed_amyloid_csv);
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
            globbed_tbl = readtable(ipr.globbed_csv, 'ReadVariableNames', true, 'Delimiter', ' '); % 3735 x 1 table
            
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
        function this = FDGDemographics(varargin)
            %% FDGDEMOGRAPHICS 
    
            this.adni_demo_ = mladni.AdniDemographics();
        end
    end
    
    %% PRIVATE
    
    properties (Access = private)
        adni_demo_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
