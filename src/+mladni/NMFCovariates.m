classdef NMFCovariates < handle
    %% Supports ADNI demographic data and other meta-data approximately at granularity of single scan objects.
    %  Prepares data for regressions more stringently than preparations for NMF
    %  
    %  Created 18-May-2022 20:19:23 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Constant)
        EXCLUSIONS = mladni.AdniDemographics.EXCLUSIONS
    end

    properties
        N_patterns
    end

    properties (Dependent)
        componentDir
        inFiles
        pet_on_T1w_suffix
        pve_1_suffix
        study_design
        targetDatasetDir
        T1w_dlicv_suffix
        weightedAverFilename
        X2mat
        X2_studydesign_mat
    end
    
    methods % GET
        function g = get.componentDir(this)
            g = fullfile(this.targetDatasetDir, sprintf('NumBases%i', this.N_patterns), 'components');
            ensuredir(g);
        end
        function g = get.inFiles(this)
            g = fullfile(this.componentDir, sprintf("%s_%s.csv", stackstr(), this.study_design));
            if ~isfile(g)
                t = table(this.table_fdg.Filelist);
                writetable(t, g, WriteVariableNames=false)
            end
        end
        function g = get.pet_on_T1w_suffix(this)
            g = this.pet_on_T1w_suffix_;
        end
        function g = get.pve_1_suffix(this)
            g = this.pve_1_suffix_;
        end
        function g = get.study_design(this)
            g = this.study_design_;
        end
        function g = get.targetDatasetDir(this)
            g = fullfile(this.data_home_, "NMF_FDG", "baseline_cn");
        end
        function g = get.T1w_dlicv_suffix(this)
            g = this.T1w_dlicv_suffix_;
        end
        function g = get.weightedAverFilename(this)
            g = fullfile(this.componentDir, ...
                sprintf("component_weighted_average_%s.csv", this.study_design));
        end
        function g = get.X2mat(this)
            g = fullfile(this.componentDir, "X2.mat");
        end
        function g = get.X2_studydesign_mat(this)
            g = fullfile(this.componentDir, sprintf("X2_%s.mat", this.study_design));
        end
    end

    methods
        function this = NMFCovariates(opts)
            %% NMFCovariates 
            % 
            % Args:
            % opts.study_design {mustBeTextScalar} = "longitudinal"
            % opts.data_home {mustBeFolder} = getenv("ADNI_HOME")
            % opts.N_patterns double = mladni.NMF.N_PATTERNS            
    
            arguments
                opts.study_design {mustBeTextScalar} = "longitudinal"
                opts.data_home {mustBeFolder} = getenv("ADNI_HOME")
                opts.N_patterns double = mladni.NMF.N_PATTERNS
            end
            this.N_patterns = opts.N_patterns;
            this.selectedNumBases_ = opts.N_patterns;
            this.study_design_ = opts.study_design;
            this.demogr_ = mladni.AdniDemographics(study_design=opts.study_design);
            this.data_home_ = opts.data_home;
            this.nmfh_ = mladni.NMFHierarchies(N_patterns = this.N_patterns);

            this.pet_on_T1w_suffix_ = 'orient-rpi_pet_on_T1w.nii.gz';
            this.T1w_dlicv_suffix_ = 'orient-rpi_T1w_dlicv.nii.gz';
            this.pve_1_suffix_ = 'orient-rpi_T1w_brain_pve_1.nii.gz';
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
            
            fdg1 = this.demogr_.table_fdg1(); % 14358 x 53 table
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
            
            fdg1 = this.demogr_.table_fdg1(); % 14358 x 53 table
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

        function f = covariates_file(this, opts)
            %% Specifies canonical f.q. filename for tables generated by NMFCovariates.
            %
            % Args:
            %
            % this mladni.NMFCovariates
            % opts.adjusted logical = false  % RStudio patterns_for_gppm.qmd has provided AdjPattern which has age, sex, apoee4 corrections
            % opts.first_scan logical = false  % by default, intended for longitudinal files
            % opts.permissive_qc logical = false  % see also this.apply_table_qc()
            % opts.selection {mustBeTextScalar} = "covariates"  % covariates, cn, preclinical, mci, ad, other
            % opts.suffix {mustBeTextScalar} = ".mat"
            %
            % Returns:
            %
            % f.q. filename

            arguments
                this mladni.NMFCovariates
                opts.adjusted logical = false  % RStudio patterns_for_gppm.qmd has provided AdjPattern which has age, sex, apoee4 corrections
                opts.first_scan logical = false  % by default, intended for longitudinal files
                opts.permissive_qc logical = false  % see also this.apply_table_qc()
                opts.selection {mustBeTextScalar} = "covariates"  % covariates, cn, preclinical, mci, ad, other
                opts.suffix {mustBeTextScalar} = ".mat"
            end

            % adjust opts.suffix
            % dbs = dbstack;
            % rname = dbs(min(2, length(dbs))).name;
            % if contains(rname, "load") || contains(rname, "save")
            %     opts.suffix = ".mat";
            % end
            % if contains(rname, "readtable") || contains(rname, "writetable")
            %     opts.suffix = ".csv";
            % end

            tag_selection = "_" + opts.selection;
            if opts.first_scan
                tag_selection = tag_selection + "_1stscan";
            end
            tag_study = "_" + this.study_design;
            fp = "NMFCovariates_table" + tag_selection + tag_study;
            if opts.adjusted
                fp = fp + "_adjpatterns";
            end
            if opts.permissive_qc
                fp = fp + "_permissive";
            end

            f = fullfile(this.componentDir, fp + opts.suffix);
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
        function t = apply_table_qc(this, t, opts)
            %% qc that is more stringent than that for NMF, needed for regressions

            arguments
                this mladni.NMFCovariates
                t table
                opts.permissive_qc logical = false
            end

            debug_file = fullfile(this.componentDir, stackstr()+this.datestr()+".mat");
            save(debug_file, 't');

            % start with 3478 rows
            
            vns = t.Properties.VariableNames;
            if any(contains(vns, 'CDGLOBAL'))
                t = t(t.CDGLOBAL ~= -1, :);
            end % 3458 rows remaining

            if opts.permissive_qc
                return
            end

            if any(contains(vns, 'Cohort'))
                t = t(t.Cohort ~= 'unknown', :);
            end % 1961 rows remaining

            if any(contains(vns, 'Components_1'))
                t = t(~isnan(t.Components_1), :);
            end % 1961 rows remaining
            if any(contains(vns, 'Age'))
                t = t(~isnan(t.Age), :);
            end % 1961 rows remaining       

%             if any(contains(vns, 'Dlicv'))
%                 t = t(t.Dlicv > 1e6, :);
%             end % 1958 rows remaining
%             if any(contains(vns, 'PVE1'))
%                 t = t(t.PVE1 > 0.3e6, :);
%             end % 1958 rows remaining
%             if any(contains(vns, 'RegErr'))
%                 t = t(t.RegErr < 2.5, :); % see also mladni.FDG()
%             end % 1941 rows remaining

            if any(contains(vns, 'ApoE4'))
                t = t(~isnan(t.ApoE4), :);
            end %

            if ~isemptytext(this.EXCLUSIONS)
                for fidx = 1:length(this.EXCLUSIONS)
                    select = ~contains(t.Filelist, this.EXCLUSIONS{fidx});
                    t = t(select, :);
                end
            end
        end
        function t = table_covariates(this, opts)
            %% Saves and returns table useful for inferences using RStudio.
            %  
            % Args:
            % this mladni.NMFCovariates
            % opts.save_1comp logical = false  % call this.table_covariates_1comp() for all components.  DEPRECATED.
            % opts.permissive_qc logical = false
            % opts.adjusted logical = false 
            % opts.reuse_cache logical = true
            %
            % N.B.:  while table_dlicv(), table_pve1(), table_regerr() may be memory-cached for all N_patterns,
            %        table_selectedComponentWeightedAverageNIFTI() varies with N_patterns, and cannot be memory-cached
            %        while supporting multiple N_patterns from a NMFCovariates object in memory.
            %        Therefor, table_covariates() should not use memory-caches, but filesystem-caches are safe so long
            %        as they are stored separately for each possible value of N_patterns.

            arguments
                this mladni.NMFCovariates
                opts.save_1comp logical = false  % call this.table_covariates_1comp() for all components.  DEPRECATED.
                opts.permissive_qc logical = false
                opts.adjusted logical = false
                opts.reuse_cache logical = true
            end

            cache_file = this.covariates_file(permissive_qc=opts.permissive_qc, adjusted=opts.adjusted);
            if opts.reuse_cache && isfile(cache_file)
                fprintf("%s: using cached from filesystem\n", stackstr())
                ld = load(cache_file);
                t = ld.t1;
                return
            end

            % from ADNIDemographics ~ 3377 rows
            t = table_fdg(this); 

            % Cohort ~ categorical
            Cohort = repmat({'unknown'}, size(t,1), 1);
            Cohort(t.CDGLOBAL == 0 & t.AmyloidStatusLong == 0) = {'CDR=0,amy-'};
            Cohort(t.CDGLOBAL == 0 & t.AmyloidStatusLong == 1) = {'CDR=0,amy+'};
            Cohort(t.CDGLOBAL == 0.5 & t.AmyloidStatusLong == 1) = {'CDR=0.5,amy+'};
            Cohort(t.CDGLOBAL > 0.5 & t.AmyloidStatusLong == 1) = {'CDR>0.5,amy+'};
            Cohort(t.CDGLOBAL > 0 & t.AmyloidStatusLong == 0) = {'CDR>0,amy-'};
            Cohort = categorical(Cohort);
            t = addvars(t, Cohort, NewVariableNames={'Cohort'});

            % DLICV ~ 3377 rows, >=4 nans
            t_ = this.table_dlicv;
            t = addvars(t, t_.Dlicv, NewVariableNames={'Dlicv'});

            % PVE1 ~ 3377 rows, >=0 nans
            t_ = this.table_pve1;
            t = addvars(t, t_.PVE1, NewVariableNames={'PVE1'});

            % RegErr ~ 3377 rows, >=29 nans
            t_ = this.table_regerr;
            t = addvars(t, t_.RegErr, NewVariableNames={'RegErr'});

            % Components ~ 3377 rows, >= 0 nans
            % implicitly marks empty Filelist for exclusion by apply_table_qc()
            t_ = this.table_selectedComponentWeightedAverageNIFTI;
            t = this.addvars_by_filelist(t, t_, t_.Components, NewVariableNames={'Components'});
            t = splitvars(t, 'Components');

            % apply table qc to obtain 3357 rows
            t = this.apply_table_qc(t, permissive_qc=opts.permissive_qc);

            % sort rows by Subject, then AcqDate
            t = sortrows(t, ["Subject", "AcqDate"]);

            % AcqDuration ~ floating-point years since first scan, after removing faulty scans
            t = this.AddAcqDuration(t);

            % save/write table
            save(cache_file, 't');
            writetable(t, strrep(cache_file, ".mat", ".csv"));

            % save separate tables for each component
            if opts.save_1comp
                for idx = 1:this.N_patterns
                    t1 = this.table_covariates_1comp(idx);
                    save(strrep(cache_file, ".mat", "_1comp.mat"), 't1');
                    writetable(t1, strrep(cache_file, ".mat", "_1comp.csv"));
                end
            end
        end
        function t = table_covariates_1comp(this, idx)
            %% returns table with indexed component in last column.

            t = this.table_covariates();
            c = t.Components;
            c = c(:, idx);
            t.Components = c;
        end
        function t = table_dlicv(this)
            %% Finds dlicv json string; estimate icv using this.find_icv().
            %  Returns:
            %      t: table with variables Filelist and DLICV.

            if ~isempty(this.table_dlicv_cache_)
                t = this.table_dlicv_cache_;
                return
            end

            Filelist = this.table_fdg.Filelist;
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
                    if ~contains(item, this.T1w_dlicv_suffix)
                        pth = fileparts(item);
                        g = glob(fullfile(pth, strcat('*', this.T1w_dlicv_suffix)));
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
            this.table_dlicv_cache_ = t;
        end
        function t = table_fdg(this)
            t = this.table_fdg5();
        end
        function t = table_fdg5(this)
            if ~isempty(this.table_fdg5_cache_)
                t = this.table_fdg5_cache_;
                return
            end
            this.table_fdg5_cache_ = this.demogr_.table_fdg5();
            t = this.table_fdg5_cache_;
        end
        function t = table_filelist(this)
            t = table(this.table_fdg.Filelist);
            t.Properties.VariableNames = {'Filelist'};
        end        
        function t = table_gppm_metarois(this)
            %% t ~ 1890x16 table

            t_ = this.table_covariates(permissive_qc=false, adjusted=true, reuse_cache=true);

            % build vars for gppm 
            RID = nan(size(t_, 1), 1);  % imputes missing MergeRid
            Time = nan(size(t_, 1), 1);  % years
            for r = 1:size(t_, 1)
                re = regexp(string(t_.Subject{r}), "\d{3}_S_(?<rid>\d{4})", "names");
                RID(r) = double(re.rid);
                Time(r) = t_.AcqDuration(r); 
            end

            % list vars from table_covariates for gppm
            vars = [ ...
                 "Metaroi", "Age", "MergeMmse", "NpBraak", "MergeCdrsb", ...
                 "MergeTau", "BRAAK1_SUVR", "BRAAK34_SUVR", "BRAAK56_SUVR", ...
                 "MergeHippocampus", "Sex", ...
                 "AmyloidStatusLong", ...
                 "SUMMARYSUVR_COMPOSITE_REFNORM", ...
                 "ApoE4", "CDGLOBAL"];
            new_vars = [ ...
                 "FDG_AD", "Age", "MMSE", "Braak_NFT", "MergeCdrsb", ...
                 "Tau_PET", "BRAAK1_SUVR", "BRAAK34_SUVR", "BRAAK56_SUVR", ...
                 "Hippo_Vol", "Sex", ...
                 "Amyloid_status", ...
                 "Amyloid_PET", ...
                 "ApoE4", "CDR"];

            % build table t
            t = table(RID, Time);
            for vidx = 1:length(vars)
                t = addvars(t, t_.(vars(vidx)), NewVariableNames=new_vars(vidx));
            end 
            t.Sex = double(strcmp(t.Sex, 'F'));

            save(this.covariates_file(selection="gppm_metarois", suffix=".mat"), "t");
            writetable(t, this.covariates_file(selection="gppm_metarois", suffix=".csv"));
        end
        function t = table_gppm_patterns(this)
            %% t ~ 1890x28 table

            t_ = this.table_covariates(permissive_qc=false, adjusted=true, reuse_cache=true);
            t_ = mergevars(t_, "AdjPattern_" + (1:24), NewVariableName="AdjPattern");  % new varname is "Var141"

            % build vars for gppm 
            RID = nan(size(t_, 1), 1);  % imputes missing MergeRid
            Time = nan(size(t_, 1), 1);  % years
            for r = 1:size(t_, 1)
                re = regexp(string(t_.Subject{r}), "\d{3}_S_(?<rid>\d{4})", "names");
                RID(r) = double(re.rid);
                Time(r) = t_.AcqDuration(r); 
            end

            % build table t
            t = table(RID, Time, t_.AdjPattern, t_.CDGLOBAL, t_.MergeCdrsb, ...
                VariableNames=["RID", "Time", "AdjPattern", "CDR", "MergeCdrsb"]);

            % split the patterns
            t = splitvars(t, "AdjPattern");
            t.Properties.VariableNames = strrep(t.Properties.VariableNames, "AdjPattern_", "Pattern");
            
            save(this.covariates_file(selection="gppm_patterns", suffix=".mat"), "t");
            writetable(t, this.covariates_file(selection="gppm_patterns", suffix=".csv"));
        end
        function t1 = table_gppm_patterns_metaroi(this)
            %% t1 ~ 1890x29 table
            % 
            %  N(CDR == 0) ~ 632
            %  N(CDR == 0.5) ~ 1052
            %  N(CDR == 1) ~ 185
            %  N(CDR == 2) ~ 18
            %  N(CDR == 3) ~ 3

            ld2 = load(this.covariates_file(selection="gppm_metarois"));
            ld1 = load(this.covariates_file(selection="gppm_patterns"));
            t1 = addvars(ld1.t, ld2.t.Amyloid_status, ld2.t.Amyloid_PET, ld2.t.FDG_AD, Before="CDR", NewVariableNames=["AmyloidStatus", "AmyloidSuvr", "Metaroi"]);
            t1.Properties.VariableNames = strrep(t1.Properties.VariableNames, "Pattern_", "Pattern");
            save(this.covariates_file(selection="gppm_patterns_metaroi", suffix=".mat"), "t1");
            writetable(t1, this.covariates_file(selection="gppm_patterns_metaroi", suffix=".csv"));
        end
        function t = table_gppm_patterns_1403(this, opts)
            %% t ~ 1403x28 table, including atypicals.  Best to remove atypicals.
            %
            % NMFCovariates_table_gppm_patterns_1403:  distribution of cdrsb:
            % cdrsb->10, n->21
            % cdrsb->9, n->3
            % cdrsb->8, n->20
            % cdrsb->7, n->24
            % cdrsb->6, n->37
            % cdrsb->5, n->62
            % cdrsb->4, n->82
            % cdrsb->3, n->118
            % cdrsb->2, n->227
            % cdrsb->1, n->413
            % cdrsb->0.5, n->286
            % cdrsb->0, n->110
            % for 872 subjects with RIDs.

            arguments
                this mladni.NMFCovariates
                opts.remove_atypical logical = true
            end

            ld = load(this.covariates_file(selection="gppm_patterns_metaroi"));
            t1 = ld.t1;

            % remove atypical dementias
            if opts.remove_atypical
                atypical = ~logical(t1.AmyloidStatus) & t1.MergeCdrsb > 0;
                t1 = t1(~atypical, :);
            end
            
            % clean variables
            t1 = removevars(t1, "AmyloidStatus");
            t1 = removevars(t1, "CDR");
            t1 = removevars(t1, "Metaroi");
            t1 = t1(~isnan(t1.MergeCdrsb), :);

            % regroup cdrsb > 7 for balanced disease phenotypes
            t1.MergeCdrsb(t1.MergeCdrsb > 7) = 7;  
            t1.MergeCdrsb(t1.MergeCdrsb == 6.5) = 6;
            t1.MergeCdrsb(t1.MergeCdrsb == 5.5) = 5;
            t1.MergeCdrsb(t1.MergeCdrsb == 4.5) = 4;
            t1.MergeCdrsb(t1.MergeCdrsb == 3.5) = 3;
            t1.MergeCdrsb(t1.MergeCdrsb == 2.5) = 2;
            t1.MergeCdrsb(t1.MergeCdrsb == 1.5) = 1;
            
            unique_cdrsb = sort(unique(t1.MergeCdrsb), 'descend');
            unique_cdrsb = unique_cdrsb(1:(end-1));

            %% find repeaters with cdrsb > 0            
            % repeaters = false(size(t1, 1), 1);
            % for cdrsb_ = asrow(unique_cdrsb)
            %     for rid_ = asrow(t1.RID(t1.MergeCdrsb == cdrsb_))
            %         repeaters = repeaters | (t1.RID == rid_);
            %     end
            % end
            % t = t1(repeaters, :);  % apply pruning filter

            t = t1;

            t.Properties.VariableNames = strrep(t.Properties.VariableNames, "Pattern", "p");
            t.Properties.VariableNames = strrep(t.Properties.VariableNames, "AmyloidSuvr", "amyloid");
            fprintf("%s:  distribution of cdrsb:\n", stackstr());
            for u_ = [asrow(unique_cdrsb), 0]
                fprintf("cdrsb->%g, n->%g\n", u_, sum(t.MergeCdrsb == u_)); 
            end
            fprintf("for %g subjects with RIDs.\n", length(unique(t.RID))); 

            Nrows = size(t, 1);
            % this.covariates_file(selection="gppm_patterns_"+Nrows, suffix=".csv")
            save(his.covariates_file(selection="gppm_patterns_"+Nrows, suffix=".mat"), "t")  
            writetable(t, his.covariates_file(selection="gppm_patterns_"+Nrows, suffix=".csv"));
        end
        function t = table_imagedataIDs(this)
            %% Finds imagedataIDs;
            %  Returns:
            %      t: table with variables Filelist and ImageDataID.

            Filelist = this.table_fdg.Filelist;
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

            if ~isempty(this.table_pve1_cache_)
                t = this.table_pve1_cache_;
                return
            end

            Filelist = this.table_fdg.Filelist;
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
                    if ~contains(item, this.pve_1_suffix)
                        pth = fileparts(item);
                        g = glob(fullfile(pth, strcat('*', this.pve_1_suffix)));
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
            this.table_pve1_cache_ = t;
        end
        function t = table_regerr(this)

            if ~isempty(this.table_regerr_cache_)
                t = this.table_regerr_cache_;
                return
            end

            Filelist = this.table_fdg.Filelist;
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
                    if ~contains(item, this.pet_on_T1w_suffix)
                        pth = fileparts(item);
                        g = glob(fullfile(pth, strcat('*', this.pet_on_T1w_suffix)));
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
            this.table_regerr_cache_ = t;
        end
        function t = table_selectedComponentWeightedAverageNIFTI(this)
            %% See also mladni.ArisCodes.calculateSelectedComponentWeightedAverageNIFTI(),
            %  which writes component_weighted_average_<study_design>.csv.
            %  t.Filelist will exclude registration failures.
            %
            %  N.B.:  caching this table in memory is incompatible with iterations of N_patterns 
            %         for a given instance of NMFCovariates.

            arguments
                this mladni.NMFCovariates
            end

            % prep param struct
            param.isList = 1 ;
            param.downSample = 1;
            param.smooth = 0;
            param.mask = fullfile(this.data_home_, "VolBin", "mask.nii.gz");
            param.permute = 0;
            param.numBase = this.N_patterns;
            param.componentDir = this.componentDir;
            assert(isfile(param.mask))

            % ensure X2.mat
            if isfile(this.X2mat)
                load(this.X2mat);
            elseif isfile(this.X2_studydesign_mat)
                load(this.X2_studydesign_mat);
            else
                data = mladni.ArisCodes.loadData(this.inFiles, param, []);
                meanX = mean(data.X, 2);
                X = data.X; % data.X(meanX > 0, :); 
                clear data;
                save(this.X2mat, "X", "meanX");
            end

            % generate table t; write table
            if ~isfile(this.weightedAverFilename)
                mladni.ArisCodes.calculateSelectedComponentWeightedAverageNIFTI( ...
                    this.inFiles, this.targetDatasetDir, this.N_patterns, this.weightedAverFilename);
            end
            t = readtable(this.weightedAverFilename, 'ReadVariableNames', false, 'Delimiter', ',');
            comp_vars = cellfun(@(x) sprintf('P%i', x), num2cell(1:this.N_patterns), UniformOutput=false);
            t.Properties.VariableNames = ['Filelist', comp_vars];
            t = mergevars(t, comp_vars, NewVariableName="Components");
        end

        %% tables of diagnostic subgroups

        function t = table_all(this, cs, T_name)
            %% N -> 1890 longitudinal, 1165 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = false
                T_name {mustBeTextScalar} = 'table_covariates'
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__;
            else
                t = this.table_firstscan(fdg__);
            end
        end
        function t = table_cn(this, varargin)
            t = this.table_cdr_0_aneg(varargin{:});
        end
        function t = table_preclinical(this, varargin)
            t = this.table_cdr_0_apos(varargin{:});
        end
        function t = table_mci(this, varargin)
            t = this.table_cdr_0p5_apos(varargin{:});
        end
        function t = table_ad(this, varargin)
            t = this.table_cdr_gt_0p5_apos(varargin{:});
        end
        function t = table_other(this, varargin)
            t = this.table_cdr_gt_0_aneg(varargin{:});
        end
        
        function t = table_cdr_0_aneg(this, cs, T_name, bl_1st)
            %% N -> 474 longitudinal, 269 cross-sectional.
            
            arguments
                this mladni.NMFCovariates
                cs logical = true
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
        function t = table_cdr_0_apos(this, cs, T_name, bl_1st)
            %% N -> 158 longitudinal, 130 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = true
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

            t1 = this.table_cdr_0_aneg(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t1.Subject)))
                t(t.Subject==dup, :) = [];
            end
        end        
        function t = table_cdr_0p5_apos(this, cs, T_name, bl_1st)
            %% N -> 533 longitudinal, 426 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = true
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

            t1 = this.table_cdr_0_aneg(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t1.Subject)))
                t(t.Subject==dup, :) = [];
            end
            t2 = this.table_cdr_0_apos(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t2.Subject)))
                t(t.Subject==dup, :) = [];
            end
        end
        function t = table_cdr_gt_0p5_apos(this, cs, T_name, bl_1st)
            %% N -> 170 longitudinal, 145 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = true
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

            t1 = this.table_cdr_0_aneg(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t1.Subject)))
                if ~isemptytext(dup)
                    t(t.Subject==dup, :) = [];
                end
            end
            t2 = this.table_cdr_0_apos(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t2.Subject)))
                if ~isemptytext(dup)
                    t(t.Subject==dup, :) = [];
                end
            end
            t3 = this.table_cdr_0p5_apos(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t3.Subject)))
                if ~isemptytext(dup)
                    t(t.Subject==dup, :) = [];
                end
            end
        end
        function t = table_cdr_gt_0_aneg(this, cs, T_name, bl_1st)
            %% N -> 555 longitudinal, 327 cross-sectional.

            arguments
                this mladni.NMFCovariates
                cs logical = true
                T_name {mustBeTextScalar} = 'table_covariates'
                bl_1st logical = false
            end

            fdg__ = this.(T_name);
            if ~cs
                t = fdg__(fdg__.CDGLOBAL > 0 & fdg__.AmyloidStatusLong == 0, :);
                return
            end
            if bl_1st
                fdg__ = this.table_firstscan(fdg__);
                t = fdg__(fdg__.CDGLOBAL > 0 & fdg__.AmyloidStatusLong == 0, :);
            else
                t = fdg__(fdg__.CDGLOBAL > 0 & fdg__.AmyloidStatusLong == 0, :);
                t = this.table_firstscan(t);
            end

            t1 = this.table_cdr_0_aneg(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t1.Subject)))
                if ~isemptytext(dup)
                    t(t.Subject==dup, :) = [];
                end
            end
            t2 = this.table_cdr_0_apos(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t2.Subject)))
                if ~isemptytext(dup)
                    t(t.Subject==dup, :) = [];
                end
            end
            t3 = this.table_cdr_0p5_apos(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t3.Subject)))
                if ~isemptytext(dup)
                    t(t.Subject==dup, :) = [];
                end
            end
            t4 = this.table_cdr_gt_0p5_apos(cs, T_name, bl_1st);
            for dup = asrow(string(intersect(t.Subject, t4.Subject)))
                if ~isemptytext(dup)
                    t(t.Subject==dup, :) = [];
                end
            end
        end 
        function t = table_firstscan(this, varargin)
            t = this.demogr_.table_firstscan(varargin{:});
        end
        
        function t = writetables(this)
            %% adds Pattern_* variables before saving

            arguments
                this mladni.NMFCovariates
            end       

            %% bootstrap:  NMFHierarchies needs NMFCovariates_table_cn_1stscan_longitudinal.mat with Components_*
            this.writetable(dx="cn", do_addvars_patterns=false);  % will be updated to include Pattern_*

            t1 = this.table_all(false);  % longitudinal
            t1 = this.nmfh_.addvars_patterns(t1);
            save(this.covariates_file(first_scan=false, suffix=".mat"), "t1");
            writetable(t1, this.covariates_file(first_scan=false, suffix=".csv"), WriteVariableNames=true);  % will overwrite

            t_cn = this.writetable(dx="cn");
            t_preclinical = this.writetable(dx="preclinical");
            t_mci = this.writetable(dx="mci");
            t_ad = this.writetable(dx="ad");
            t_other = this.writetable(dx="other");

            t = [t_cn; t_preclinical; t_mci; t_ad; t_other];
            t_1st = this.table_all(true);  % simple selection of 1st scan of FDG
            assert(size(t_1st, 1) == ...
                size(t_cn, 1) + size(t_preclinical, 1) + size(t_mci, 1) + size(t_ad, 1) + size(t_other, 1))
            t = this.nmfh_.addvars_patterns(t);
            save(this.covariates_file(first_scan=true, suffix=".mat"), "t");
            writetable(t, this.covariates_file(first_scan=true, suffix=".csv"), WriteVariableNames=true);  % will overwrite with corrections
        end
        function t = writetable(this, opts)
            %% adds Pattern_* variables before saving
            
            arguments
                this mladni.NMFCovariates
                opts.dx {mustBeTextScalar} = "cn"
                opts.first_scan logical = true
                opts.do_addvars_patterns logical = true;  
            end

            if contains(opts.dx, "cn", IgnoreCase=true) || contains(opts.dx, "cdr=0,amy-", IgnoreCase=true)
                t = this.table_cn(opts.first_scan);
            end
            if contains(opts.dx, "preclinical", IgnoreCase=true) || contains(opts.dx, "cdr=0,amy+", IgnoreCase=true)
                t = this.table_preclinical(opts.first_scan);
            end
            if contains(opts.dx, "mci", IgnoreCase=true) || contains(opts.dx, "cdr=0.5,amy+", IgnoreCase=true)
                t = this.table_mci(opts.first_scan);
            end
            if contains(opts.dx, "ad", IgnoreCase=true) || contains(opts.dx, "cdr>0.5,amy+", IgnoreCase=true)
                t = this.table_ad(opts.first_scan);
            end
            if contains(opts.dx, "other", IgnoreCase=true) || contains(opts.dx, "cdr>0,amy-", IgnoreCase=true)
                t = this.table_other(opts.first_scan);
            end
            
            if opts.do_addvars_patterns
                t = this.nmfh_.addvars_patterns(t);
            end
            save(this.covariates_file(selection=opts.dx, first_scan=true, suffix=".mat"), "t");
            writetable(t, this.covariates_file(selection=opts.dx, first_scan=true, suffix=".csv"), WriteVariableNames=true)

            fprintf("Age range (years at enrollment):  %g - %g\n", min(t.Age), max(t.Age));
            fprintf("Age (years at enrollment):  mu = %g, sigma = %g\n", mean(t.Age), std(t.Age));
            fprintf("Female (pcnt):  %g\n", sum(strcmp(t.Sex,"F"))/size(t,1));
            fprintf("No. APOEe4 range (alleles):  %g - %g\n", min(t.ApoE4), max(t.ApoE4));
            fprintf("No. APOEe4 (alleles):  mu = %g, sigma = %g\n", mean(t.ApoE4), std(t.ApoE4));
        end
    end

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
                    if ~contains(item, this.T1w_dlicv_suffix)
                        pth = fileparts(item);
                        g = glob(fullfile(pth, strcat('*', this.T1w_dlicv_suffix)));
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
                    if ~contains(item, this.pve_1_suffix)
                        pth = fileparts(item);
                        g = glob(fullfile(pth, strcat('*', this.pve_1_suffix)));
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
        function str = datestr()
            str = string(datetime("now", Format="yyyyMMddhhmmss"));
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
    
    %% PROTECTED
    
    properties (Access = protected)
        data_home_
        demogr_
        nmfh_
        pet_on_T1w_suffix_
        pve_1_suffix_
        selectedNumBases_
        study_design_
        table_cohorts_
        table_dlicv_cache_
        table_fdg5_cache_
        table_pve1_cache_
        table_regerr_cache_
        T1w_dlicv_suffix_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
