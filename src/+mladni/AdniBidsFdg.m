classdef AdniBidsFdg < mladni.AdniBids
    %% 
    %  
    %  Created 10-Feb-2022 19:47:37 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1837725 (R2021b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function [j,c] = parcluster()
            %% #PARCLUSTER2
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/          
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.AdniBidsFdg.batch, 1, {'len', []}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false); % {'len', []}, 'Pool', 31            
        end
        function t = batch(varargin)
            %% From subjectsDir, dcm -> nifti or nifti -> nifti, for ADNI preprocessed
            %  Params:
            %      len (numeric):  # subjects
            %      subjectsDir (folder):  loni folders extracted from archives
            %      bidsDir (folder):  default $SINGULARITY_HOME?ADNI/bids/rawdata
            
            disp('Start mladni.AdniBidsFdg.batch()')            
            t0 = tic;

            subjectsDir = fullfile('/home/aris_data', 'ADNI_FDG', 'loni_FDG', 'loni', '');
            bidsDir = fullfile('/home/aris_data', 'ADNI_FDG', 'bids', 'rawdata', '');
            
            ip = inputParser;
            addParameter(ip, 'len', [], @isnumeric)
            addParameter(ip, 'subjectsDir', subjectsDir, @isfolder)
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:});
            ipr = ip.Results;      
            
            pwd0 = pushd(ipr.subjectsDir);
            globbed = globFoldersT(fullfile(ipr.subjectsDir, '*_S_*'));
            if ~isempty(ipr.len)
                globbed = globbed(1:ipr.len);
            end
            fprintf('mladni.AdniBidsFdg.batch.globbed:\n')
            disp(size(globbed))
            
            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    im = globFoldersT(fullfile(globbed{idx}, 'Coreg*/*/I*'));                    
                    this = mladni.AdniBidsFdg();
                    for i = 1:length(im)
                        this.build_nifti(im{i}, ipr.bidsDir); %#ok<PFBNS>
                        this.update_json(im{i}, ipr.bidsDir);
                    end
                catch ME
                    handwarning(ME)
                end
            end
            popd(pwd0);
            
            t = toc(t0);
            disp('mladni.AdniBidsT1w.batch() completed')  
        end
        function create(varargin)
            %% creates, e.g., /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-20090713133153_trc-FDG_proc-CAS_pet.json
            %                 /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-20090713133153_trc-FDG_proc-CAS_pet.nii.gz
            %  Args:
            %      subjectsDir (folder): e.g., $SINGULARITY_HOME/ADNI/loni_FDG/loni; default pwd.
            %      bidsDir (folder): default $SINGULARITY_HOME?ADNI/bids/rawdata

            ip = inputParser;
            subjectsDir = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'loni_FDG', 'loni', '');
            bidsDir = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'rawdata', '');
            addParameter(ip, 'subjectsDir', subjectsDir, @isfolder)
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            pwd0 = pushd(ipr.subjectsDir);
            this = mladni.AdniBidsFdg();
            s = globFoldersT(fullfile(ipr.subjectsDir, '*_S_*'));
            for si = 1:length(s)
                im = globFoldersT(fullfile(s{si}, 'Coreg*/*/I*'));
                for i = 1:length(im)
                    this.build_nifti(im{i}, ipr.bidsDir);
                    this.update_json(im{i}, ipr.bidsDir);
                end
            end
            popd(pwd0);
        end
        function [tbl,csv_fn] = create_table_fdg_t1w(varargin)
            % Args:
            %      bidsDir (folder): default $SINGULARITY_HOME?ADNI/bids/rawdata
            % Returns:
            %      tbl: table for Matlab pairing fdg with best t1w for warping to atlases.
            %      csv_fn: exported for use in other platforms.

            ip = inputParser;
            bidsDir = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'rawdata', '');
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            pwd0 = pushd(ipr.bidsDir);
            this = mladni.AdniBidsFdg();
            fdgs = [];
            t1ws = [];

            s = globFoldersT(fullfile(ipr.bidsDir, 'sub-*'));
            for si = 1:length(s)

                % for each subject
                fdgs_ = glob(fullfile(s{si}, 'ses-*', 'pet', '*.nii.gz'));
                t1ws_ = glob(fullfile(s{si}, 'ses-*', 'anat', '*.nii.gz'));
                t1ws__ = cell(size(fdgs_));
                for f = 1:length(fdgs_)
                    t1ws__{f} = this.select_t1w(t1ws_, fdgs_{f});  
                end

                fdgs = vertcat(fdgs, fdgs_); %#ok<AGROW> 
                t1ws = vertcat(t1ws, t1ws__); %#ok<AGROW> 
            end
            popd(pwd0)
            
            % assemble table
            tbl = table(fdgs, t1ws, 'VariableNames', {'FDG_filename', 'T1w_filename'});
            csv_fn = fullfile(ipr.bidsDir, 'rosters', mladni.AdniBidsFdg.TABLE_FDG_T1W_FILE);
            writetable(tbl, csv_fn);
        end
    end

    properties (Constant)
        TABLE_FDG_T1W_FILE = 'mladni_AdniBidsFdg_table_fdg_t1w.csv'
    end

    methods
        function fn = build_nifti(~, image_folder, rawdata_folder)
            %% calls dcm2niix
            %  Args:
            %      image_folder (folder): e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0/I12345/
            %      rawdata_folder (folder): e.g., /path/to/bids/rawdata
            %  Returns
            %      fn: e.g., /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-yyyyMMddHHmmss_trc-fdg_proc-CAS_pet.nii.gz

            assert(strcmp(image_folder(1), filesep)) % abs path

            % sub
            dt_folder = fileparts(image_folder); % e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0
            [~,sub_str] = fileparts(fileparts(fileparts(dt_folder)));
            sub = ['sub-' strrep(sub_str, '_', '')];

            % ses
            [~,dt_str] = fileparts(dt_folder); % e.g., yyyy-MM-dd_HH_mm_ss.0
            dt = datetime(strtok(dt_str, '.'), 'InputFormat', 'yyyy-MM-dd_HH_mm_ss');
            ses8 = sprintf('ses-%i%02i%02i', ...
                year(dt), month(dt), day(dt));
            ses14 = sprintf('ses-%i%02i%02i%02i%02i%02i', ...
                year(dt), month(dt), day(dt), hour(dt), minute(dt), second(dt));

            % proc
            [~,proc_str] = fileparts(fileparts(dt_folder));
            proc = '';
            if startsWith(proc_str, 'co', 'IgnoreCase', true)
                proc = 'C';
            end
            if contains(proc_str, 'dynamic', 'IgnoreCase', true)
                proc = [proc 'D'];
            end
            if contains(proc_str, 'av', 'IgnoreCase', true)
                proc = [proc 'A'];
            end
            if contains(proc_str, 'vox', 'IgnoreCase', true)
                proc = [proc 'S'];
            end
            if contains(proc_str, 'uniform', 'IgnoreCase', true)
                proc = [proc 'U'];
            end
            if isempty(proc)
                proc = 'unknown';
            end
            proc = ['proc-' proc];

            % dcm2niix
            dest_folder = fullfile(rawdata_folder, sub, ses8, 'pet', '');
            if ~isfolder(dest_folder)
                mkdir(dest_folder);
            end
            base = sprintf('%s_%s_trc-FDG_%s_pet', sub, ses14, proc);            
            if ~isfile(fullfile(dest_folder, [base '.nii.gz']))
                mlpipeline.Bids.dcm2niix(image_folder, 'f', base, 'o', dest_folder)
            end
            fn = fullfile(dest_folder, strcat(base, '.nii.gz'));
        end
        function p = quality(this, rid, dt)
            t = this.table_qc_(this.table_qc_.RID == rid & this.table_qc_.EXAMDATE == dt, :);
            if ~isempty(t)
                p = t.PASS;
                return
            end
            t = this.table_c3_(this.table_c3_.RID == rid & this.table_c3_.SCANDATE == dt, :);
            if ~isempty(t)
                p = t.SCANQLTY;
                return
            end
            p = nan;
        end
        function p = phase(this, rid, dt)
            t = this.table_qc_(this.table_qc_.RID == rid & this.table_qc_.EXAMDATE == dt, :);
            if ~isempty(t)
                p = t.Phase;
                return
            end
            t = this.table_c3_(this.table_c3_.RID == rid & this.table_c3_.SCANDATE == dt, :);
            if ~isempty(t)
                p = t.Phase;
                return
            end
            p = 'UNKNOWN';
        end
        function fn = select_t1w(this, t1w_list, fdg_fn)
            %  Args:
            %      t1w_list (cell): fqfilenames for a subject.
            %      fdg_fn(file): fqfilename for a subject.
            %  Returns:
            %      fn: selected from t1w_list.

            % prune t1w_list for more processing
            t1w_list = t1w_list(~contains(t1w_list, '-mask'));
            t1w_list = this.select_t1w_quality(t1w_list);

            % select closest 
            try
                t1w_dt = this.fqnifti2dt(t1w_list);
                fdg_dt = this.nii2datetime(basename(fdg_fn));
                min_ = min(abs(t1w_dt - fdg_dt));
                t1w_list = t1w_list(abs(t1w_dt - fdg_dt) == min_);

                lens_ = cell2mat(cellfun(@(x) length(x), t1w_list, 'UniformOutput', false)); % lengths of proc identifying strings
                [~,idx] = max(lens_); % longer filename has more processing
                fn = t1w_list{idx};
            catch ME
                handwarning(ME)
                fn = t1w_list{1}; % earliest available t1w, which should have least pathology
            end
        end
        function t1w_list = select_t1w_quality(this, t1w_list)
            %  Args:
            %      t1w_list (cell): t1w fqfilenames for a subject

            rid = this.sub2rid(this.nifti2sub(mybasename(t1w_list{1})));
            q = this.table_mri_quality_(this.table_mri_quality_.RID == rid, :);
            if all(q.PASS == 1)
                return
            end

            q = q(q.PASS ~= 1, :);
            defects = false(size(t1w_list));
            for qi = 1:size(q,1)
                defects = defects | ...
                    (this.fqnifti2dt(t1w_list) == q.EXAMDATE(qi) & ...
                     strcmp(this.fqnifti2uid(t1w_list), q.LONIUID(qi)));
            end
            t1w_list = t1w_list(~defects);
        end
        function fn = update_json(~, image_folder, rawdata_folder)
            %% calls jsonrecode
            %  Args:
            %      image_folder (folder): e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0/I12345/
            %      rawdata_folder (folder): e.g., /path/to/bids/rawdata
            %  Returns
            %      fn: e.g., /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-yyyyMMddHHmmss_trc-fdg_proc-CAS_pet.nii.gz

            assert(strcmp(image_folder(1), filesep)) % abs path

            % sub
            dt_folder = fileparts(image_folder); % e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0
            [~,sub_str] = fileparts(fileparts(fileparts(dt_folder)));
            sub = ['sub-' strrep(sub_str, '_', '')];

            % ses
            [~,dt_str] = fileparts(dt_folder); % e.g., yyyy-MM-dd_HH_mm_ss.0
            dt = datetime(strtok(dt_str, '.'), 'InputFormat', 'yyyy-MM-dd_HH_mm_ss');
            ses8 = sprintf('ses-%i%02i%02i', ...
                year(dt), month(dt), day(dt));
            ses14 = sprintf('ses-%i%02i%02i%02i%02i%02i', ...
                year(dt), month(dt), day(dt), hour(dt), minute(dt), second(dt));

            % proc
            [~,proc_str] = fileparts(fileparts(dt_folder));
            proc = '';
            if startsWith(proc_str, 'co', 'IgnoreCase', true)
                proc = 'C';
            end
            if contains(proc_str, 'dynamic', 'IgnoreCase', true)
                proc = [proc 'D'];
            end
            if contains(proc_str, 'av', 'IgnoreCase', true)
                proc = [proc 'A'];
            end
            if contains(proc_str, 'vox', 'IgnoreCase', true)
                proc = [proc 'S'];
            end
            if contains(proc_str, 'uniform', 'IgnoreCase', true)
                proc = [proc 'U'];
            end
            if isempty(proc)
                proc = 'unknown';
            end
            proc = ['proc-' proc];

            % jsonrecode
            dest_folder = fullfile(rawdata_folder, sub, ses8, 'pet');
            base = sprintf('%s_%s_trc-FDG_%s_pet', sub, ses14, proc);
            j = fullfile(dest_folder, strcat(base, '.json'));
            dcms = glob(fullfile(image_folder, '*.dcm'));
            s.ADNI_INFO.ScanDate = datestr(dt, 'yyyy-mm-dd');
            s.ADNI_INFO.OriginalPath = image_folder;
            s.ADNI_INFO.OriginalFilename = basename(dcms{1});
            jsonrecode(j, s, 'filenameNew', j);
        end

        function this = AdniBidsFdg(varargin)
            this = this@mladni.AdniBids(varargin{:});
            
            this.ad_ = mladni.AdniDemographics();      
            this.table_mri_quality_ = this.ad_.table_mri_quality();
            %this.t1_filenames_ = readtable(fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'rawdata', 'rosters', 't1_filenames.csv'));
            %assert(~isempty(this.t1_filenames_))
            %this.table_qc_ = this.ad_.table_pet_qc();
            %this.table_c3_ = this.ad_.table_pet_c3();
        end
    end
    
    %% PROTECTED

    properties (Access = protected)
        ad_
        table_mri_quality_
        table_qc_
        table_c3_
        t1_filenames_
    end

    methods (Access = protected)
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
