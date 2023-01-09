classdef FDGQC < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% line1
    %  
    %  num_scans(conditions) ~ [1025, 1647, 794], sum ~ 3466
    %  num_scans regardless of known conditions ~ 3735
    %  num_scans(completed_fdg_warped) ~ 3709
    %  num_scans(t4_resolve_rerr <= 10) ~ 3560
    %  num_scans(t4_resolve_terr <= 10) ~ 3566
    %  
    %  Created 21-Mar-2022 14:24:57 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    methods (Static)
        function propcluster()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '10000'; % deepmrseg requires 10 GB; else 5 GB
            c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = '';
            c.AdditionalProperties.WallTime = '24:00:00'; % 24 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function getDebugLog(j,c)
            try
                c.getDebugLog(j)
            catch
                c.getDebugLog(j.Tasks(end))
            end
        end
        function [j,c] = parcluster()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.FDGQC.batch, 1, {}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);            
        end
        function [j,c] = parcluster2()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.FDGQC.batch2, 1, {}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);            
        end        
        function [j,c] = par_build_geom_stats()
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.FDGQC.batch_build_geom_stats, 1, {}, 'Pool', 9, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);
        end
        function [j,c] = par_build_fast_stats()
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.FDGQC.batch_build_fast_stats, 1, {}, 'Pool', 9, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);
        end
        function [j,c] = par_find_incomplete()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.FDGQC.batch_find_incomplete, 1, {}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);            
        end        
        
        function t = batch(varargin)
            %% for all globbed:  
            %      find t4_resolve cost_final as rerr, terr
            %  save rerr.mat, terr.mat, fcost.mat
            
            ip = inputParser;
            addParameter(ip, 'proc', mladni.FDG.PROC, @istext)
            addParameter(ip, 'tag', '_orient-rpi', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mladni.FDG.batch()')        
            
            t0 = tic;            
            setenv('ADNI_HOME', '/scratch/jjlee/Singularity/ADNI')
            globbed = globT( ...
                fullfile( ...
                    getenv('ADNI_HOME'), 'bids', 'derivatives', 'sub-*', 'ses-*', 'pet', ...
                    sprintf('*-%s%s_pet_on_T1w.json', ipr.proc, ipr.tag)));
            fprintf('mladni.FDG.batch.globbed.size:\n')
            disp(size(globbed))
            rerr = nan(size(globbed));
            terr = nan(size(globbed));
            
            parfor idx = 1:length(globbed)
                try
                    j = jsondecode(fileread(globbed{idx}));
                    rerr(idx) = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_rotation_error.err); %#ok<PFOUS>
                    terr(idx) = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_translation_error.err); %#ok<PFOUS>
                catch ME
                    handwarning(ME);
                end
            end
            save(fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', 'rerr.mat'), 'rerr')
            save(fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', 'terr.mat'), 'terr')
            t = toc(t0);

            disp('mladni.FDG.batch() completed')            
        end
        function t = batch2(varargin)
            %% for all FDG with resolve error > err:
            %      symlink .nii.gz, .json to derivatives/QC/resolve_error_gt_error
            %
            %  Args:
            %      err (scalar):  default is 20 degrees/mm.
            
            ip = inputParser;
            addParameter(ip, 'proc', mladni.FDG.PROC, @istext)
            addParameter(ip, 'tag', '_orient-rpi-std', @istext)
            addParameter(ip, 'err', 20, @isscalar)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mladni.FDG.batch()')        
            
            t0 = tic;            
            setenv('ADNI_HOME', '/scratch/jjlee/Singularity/ADNI')
            qc_path = sprintf('%s/bids/derivatives/QC/resolve_error_gt_%g', getenv('ADNI_HOME'), ipr.err);
            ensuredir(qc_path)
            globbed = globT( ...
                fullfile( ...
                    getenv('ADNI_HOME'), ...
                    'bids/derivatives', ...
                    'sub-*', 'ses-*', 'pet', ...
                    sprintf('*-%s%s_pet_final.json', ipr.proc, ipr.tag)));
            fprintf('mladni.FDG.batch.globbed.size:\n')
            disp(size(globbed))
            
            parfor idx = 1:length(globbed)
                try
                    j = jsondecode(fileread(globbed{idx}));
                    rerr = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_rotation_error.err);
                    terr = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_translation_error.err);
                    if rerr > ipr.err && terr > ipr.err %#ok<PFBNS>
                        json = globbed{idx};
                        niigz = strrep(json, '.json', '.nii.gz')
                        mlbash(sprintf('ln -s %s %s', json, qc_path));
                        mlbash(sprintf('ln -s %s %s', niigz, qc_path));
                    end
                catch ME
                    handwarning(ME);
                end
            end
            t = toc(t0);

            disp('mladni.FDG.batch() completed')            
        end
        function t = batch_find_incomplete(varargin)
            
            t0 = tic;
            derpth = fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', '');
            subpth = glob(fullfile(derpth, 'sub-*', ''));
            len = length(subpth);
            gs = cell(len, 1);
            parfor idx = 1:len                
                % each subject path
                
                petpth = glob(fullfile(subpth{idx}, 'ses-*', 'pet', ''));
                c = {};
                for ip = 1:length(petpth)
                    g = glob(fullfile(petpth{ip}, ...
                        'sub-*_ses-*_trc-FDG_proc-CASU_orient-rpi-std_pet_final.nii.gz'));
                    if isempty(g)
                        c = vertcat(c, petpth{ip});
                    end
                end
                gs{idx} = c;
            end
            save(fullfile(derpth, 'gs.mat'), 'gs');
            gs = gs(cellfun(@(x) ~isempty(x), gs));
            gs = vertcat(gs{:});
            tbl = cell2table(gs, 'VariableNames', {'pet_path_incomplete'});
            writetable(tbl, fullfile(derpth, 'pet_paths_incomplete.csv'));            
            t = toc(t0);
            
            disp('mladni.FDG.batch_find_incomplete() completed')
        end 
        function t = batch_build_geom_stats(varargin)
            %% 
            
            ip = inputParser;
            addParameter(ip, 'new_tag', '_on_T1w_Warped', @istext); % '_on_T1w_Warped'
            addParameter(ip, 'proc', strcat('proc-',mladni.FDG.PROC,'-ponsvermis'), @istext);
            addParameter(ip, 'conditions', {'cn', 'mci', 'dementia', ...
                                            'cn_amypos', 'mci_amypos', 'dementia_amypos', ...
                                            'cn_amyneg', 'mci_amyneg', 'dementia_amyneg'}, @iscell);
            parse(ip, varargin{:});
            ipr = ip.Results;
            new_tag = ipr.new_tag;
            proc = ipr.proc;
            conditions = ipr.conditions;
            
            disp('Start mladni.FDGQC.batch_build_geom_stats()')  
            
            t0 = tic;
            c = cell(1, length(conditions));
            fileprefix_mean = strcat('sub-all_ses-all_trc-FDG_', proc, '_orient-rpi_pet', new_tag, '_mean');
            fileprefix_var = strcat('sub-all_ses-all_trc-FDG_', proc, '_orient-rpi_pet', new_tag, '_var');
            fileprefix_std = strcat('sub-all_ses-all_trc-FDG_', proc, '_orient-rpi_pet', new_tag, '_std');
            
            parfor idx = 1:length(conditions)
                try         
                    mladni.CHPC3.setenvs();
                    csv_file = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', ...
                        strcat('globbed_', conditions{idx},'.csv'));
                    fprintf('########## mladni.FDGQC.batch_build_geom_stats().csv_file -> %s ##########\n', csv_file)
                    assert(isfile(csv_file));
                    tbl = readtable(csv_file, 'ReadVariableNames', false, 'Delimiter', ' ');
                    fns = strrep(tbl.(['rawdata_pet_filename_' conditions{idx}]), ...
                        'rawdata', 'derivatives');
                    fns = strrep(fns, ...
                        '.nii.gz', strcat(new_tag, '.nii.gz'));
                    fns = strrep(fns, ...
                        strcat('proc-', mladni.FDG.PROC), proc);
                    len = length(fns);
                    
                    % mean
                    ic = mlfourd.ImagingContext2(fns{1});
                    fprintf('read %s\n', fns{1});
                    e_count = 0;
                    for ifn = 2:len
                        try
                            ic_ = mlfourd.ImagingContext2(fns{ifn});
                            assert(dipsum(isnan(ic_)) == 0)
                            ic = ic + ic_;
                            ic.fileprefix = fileprefix_mean;
                            %fprintf('read %s\n', fns{ifn});
                        catch ME
                            e_count = e_count + 1;
                            handwarning(ME);
                        end
                    end
                    len = len - e_count;
                    ic = ic/len;
                    ic.filepath = myfileparts(csv_file);
                    ic.fileprefix = fileprefix_mean;
                    c{idx}.mean = ic;
                    c{idx}.mean_len = len;
                    c{idx}.mean_exception_count = e_count;
                    
                    % var
                    ic1 = (mlfourd.ImagingContext2(fns{1}) - c{idx}.mean).^2;
                    for ifn = 2:len
                        try
                            ic1_ = mlfourd.ImagingContext2(fns{ifn});
                            assert(dipsum(isnan(ic1_)) == 0)
                            ic1 = ic1 + (ic1_ - c{idx}.mean).^2;
                            ic1.fileprefix = fileprefix_var;                        
                        catch ME
                            e_count = e_count + 1;
                            handwarning(ME);
                        end
                    end
                    len = len - e_count;
                    ic1 = ic1/(len - 1);
                    ic1.filepath = myfileparts(csv_file);
                    ic1.fileprefix = fileprefix_var;
                    c{idx}.var = ic1;
                    c{idx}.var_len = len;
                    c{idx}.var_exception_count = e_count;
                    
                    % std
                    c{idx}.std = sqrt(c{idx}.var);
                    c{idx}.std.filepath = myfileparts(csv_file);
                    c{idx}.std.fileprefix = fileprefix_std;
                catch ME
                    handwarning(ME);
                end
            end   
            
            save( ...
                fullfile('/scratch', 'jjlee', 'Singularity', ...
                'ADNI', 'bids', 'derivatives', strcat('FDGQC_batch_build_geom_stats', new_tag, '_', proc, '.mat')), 'c');
            t = toc(t0);
            
            disp('mladni.FDGQC.batch_build_geom_stats() completed')   
        end
        function t = batch_build_fast_stats(varargin)
            %% 
            
            disp('Start mladni.FDGQC.batch_build_fast_stats()')  
            
            t0 = tic;            
            conditions = {'pve0_cn', 'pve0_mci', 'pve0_dementia', ...
                          'pve1_cn', 'pve1_mci', 'pve1_dementia', ...
                          'pve2_cn', 'pve2_mci', 'pve2_dementia'};
            c = cell(1, 9);
            new_tag = '_detJ_Warped'; 
            
            parfor idx = 1:length(conditions)
                try         
                    mladni.CHPC3.setenvs();
                    
                    fileprefix_mean = strcat('sub-all_ses-all_orient-rpi_T1w_brain_', conditions{idx}, '_mean');
                    fileprefix_var = strcat('sub-all_ses-all_orient-rpi_T1w_brain_', conditions{idx}, '_var');
                    fileprefix_std = strcat('sub-all_ses-all_orient-rpi_T1w_brain_', conditions{idx}, '_std');
                    csv_file = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', ...
                        strcat('globbed_', conditions{idx},'.csv'));
                    fprintf('########## mladni.FDGQC.batch_build_fast_stats().csv_file -> %s ##########\n', csv_file)
                    assert(isfile(csv_file));
                    tbl = readtable(csv_file, 'ReadVariableNames', false, 'Delimiter', ' ');
                    fns = tbl.(['filename_' conditions{idx}]);
                    len = length(fns);
                    
                    % mean
                    ic = mlfourd.ImagingContext2(fns{1});
                    fprintf('read %s\n', fns{1});
                    e_count = 0;
                    for ifn = 2:len
                        try
                            ic = ic + mlfourd.ImagingContext2(fns{ifn});
                            ic.fileprefix = fileprefix_mean;
                            %fprintf('read %s\n', fns{ifn});
                        catch ME
                            e_count = e_count + 1;
                            handwarning(ME);
                        end
                    end
                    len = len - e_count;
                    ic = ic/(len - 1);
                    ic.filepath = myfileparts(csv_file);
                    ic.fileprefix = fileprefix_mean;
                    c{idx}.mean = ic;
                    c{idx}.mean_len = len;
                    c{idx}.mean_exception_count = e_count;
                    
                    % var
                    ic1 = (mlfourd.ImagingContext2(fns{1}) - c{idx}.mean).^2;
                    for ifn = 2:len
                        try
                            ic1 = ic1 + (mlfourd.ImagingContext2(fns{ifn}) - c{idx}.mean).^2;
                            ic1.fileprefix = fileprefix_var;
                        catch ME
                            e_count = e_count + 1;
                            handwarning(ME);
                        end
                    end
                    len = len - e_count;
                    ic1 = ic1/(len - 1);
                    ic1.filepath = myfileparts(csv_file);
                    ic1.fileprefix = fileprefix_var;
                    c{idx}.var = ic1;
                    c{idx}.var_len = len;
                    c{idx}.var_exception_count = e_count;
                     
                    % std
                    c{idx}.std = sqrt(c{idx}.var);
                    c{idx}.std.filepath = myfileparts(csv_file);
                    c{idx}.std.fileprefix = fileprefix_std;
                catch ME
                    handwarning(ME);
                end
            end   
            
            save( ...
                fullfile('/scratch', 'jjlee', 'Singularity', ...
                'ADNI', 'bids', 'derivatives', strcat('batch_build_fast_stats', new_tag, '.mat')), 'c');
            t = toc(t0);
            
            disp('mladni.FDGQC.batch_build_fast_stats() completed')   
        end
        function check_conservation_of_fdg_activity()
        end
        function check_conservation_of_dlicv_mass()
        end
        function check_conservation_of_fast_mass()
        end        
        function check_positivity_of_fdg_warped()
        end
        function t = foo(varargin)
            t0 = tic;
            trash = '/scratch/jjlee/Singularity/ADNI/.Trash';
            ensuredir(trash);
            csv = fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', 'pet_paths_incomplete.csv');
            tbl = readtable(csv, 'Delimiter', ',');
            for it = 1:size(tbl,1)
                pth = tbl.pet_path_incomplete{it};
                src = strrep(pth, getenv('ADNI_HOME'), '/scratch/jjlee/Singularity/ADNI');
                dest = strrep(pth, getenv('ADNI_HOME'), '/scratch/jjlee/Singularity/ADNI/.Trash');
                movefile(src, dest);
            end
            t = toc(t0);
        end
        function histograms_final_costs()
            f = load('fcost.mat'); % 12-affine corratio, T1w \rightarrow atlas
            figure;
            histogram(f.fcost, "NumBins", 100);
            xlabel("12-affine corratio, T1w \rightarrow atlas");
            saveFigures(pwd, 'closeFigure', true, 'prefix', 'fcost');
            
            r = load('rerr.mat'); % rotational err, degrees, T1w \rightarrow atlas
            figure;
            histogram(r.rerr, "NumBins", 100);
            xlabel("rotational err, degrees, PET \rightarrow T1w");
            saveFigures(pwd, 'closeFigure', true, 'prefix', 'rerr');
            
            t = load('terr.mat'); % translational err, mm, T1w \rightarrow atlas
            figure;
            histogram(t.terr, "NumBins", 100);
            xlabel("translational err, mm, PET \rightarrow T1w");
            saveFigures(pwd, 'closeFigure', true, 'prefix', 'terr');
        end
    end
    
    properties (Dependent)
        filepath
    end
    
    methods
        
        %% GET
        
        function g = get.filepath(~)
            g = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', '');
        end
        
        %%
        
        function this = FDGQC(varargin)
            %%
            %  Args:
            %      load_tables (logical)
            
            ip = inputParser;
            addParameter(ip, "load_tables", true, @islogical);
            parse(ip, varargin{:})
            ipr = ip.Results;  
            
            try
                if ipr.load_tables
                    this.table_cn();
                    this.table_mci();
                    this.table_dementia();
                end            
            catch ME
                handwarning(ME)
            end
        end
        
        function t = table_cn(this)
            if ~isempty(this.table_cn_)
                t = this.table_cn_;
                return
            end
            fn = fullfile(this.filepath, 'globbed_cn.csv');
            this.table_cn_ = readtable(fn, 'ReadVariableNames', false, 'Delimiter', ' ');
            t = this.table_cn_;
        end
        function t = table_mci(this)
            if ~isempty(this.table_mci_)
                t = this.table_mci_;
                return
            end
            fn = fullfile(this.filepath, 'globbed_mci.csv');
            this.table_mci_ = readtable(fn, 'ReadVariableNames', false, 'Delimiter', ' ');
            t = this.table_mci_;
        end
        function t = table_dementia(this)
            if ~isempty(this.table_dementia_)
                t = this.table_dementia_;
                return
            end
            fn = fullfile(this.filepath, 'globbed_dementia.csv');
            this.table_dementia_ = readtable(fn, 'ReadVariableNames', false, 'Delimiter', ' ');
            t = this.table_dementia_;
        end
        
        function t = table_dx_amyloid(this, varargin)
            ip = inputParser;
            addRequired(ip, 'dx', @istext);
            addOptional(ip, 'amyloid_status', nan, @isscalar)
            parse(ip, varargin{:});        
            ipr = ip.Results;
            ipr.dx = lower(ipr.dx);
            
            switch ipr.dx
                case 'cn'
                    t = this.table_cn();
                case 'mci'
                    t = this.table_mci();
                case 'dementia'
                    t = this.table_dementia();
                otherwise
                    error('mladni:RuntimeError', 'FDGQC.table_dx_amyloid');
            end
            if isnan(ipr.amyloid_status) || ipr.amyloid_status < 0 
                return
            end
            
            ad = mladni.AdniDemographics();
            tad = ad.table_fdg1();
            as = tad.AmyloidStatus;
            as(isnan(as)) = -1;
            sub = strrep(tad.Subject, '_', '');
            dt = string(datestr(tad.AcqDate, 'yyyymmdd'));
            if 0 == ipr.amyloid_status
                v = strcat(t.Properties.VariableNames{1}, '_amyneg');
                cells = t{:,:};
                rid_ = sub(0 == as);
                dt_ = dt(0 == as);
                select = cellfun(@(x) contains(x, rid_) && contains(x, dt_), cells);
                t = table(cells(select));
                t.Properties.VariableNames = {v};
                writetable(t, fullfile(this.filepath, strcat('globbed_', ipr.dx, '_amyneg.csv')));
            end
            if ipr.amyloid_status > 0
                v = strcat(t.Properties.VariableNames{1}, '_amypos');
                cells = t{:,:};
                rid_ = sub(1 == as);
                dt_ = dt(1 == as);
                select = cellfun(@(x) contains(x, rid_) && contains(x, dt_), cells);
                t = table(cells(select));
                t.Properties.VariableNames = {v};
                writetable(t, fullfile(this.filepath, strcat('globbed_', ipr.dx, '_amypos.csv')));
            end
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)  
        table_cn_
        table_mci_
        table_dementia_
    end    
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
