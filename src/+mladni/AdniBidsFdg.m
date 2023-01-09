classdef AdniBidsFdg < handle & mladni.AdniBidsT1w
    %% AdniBidsFdg migrates imaging from loni_FDG folders to bids/rawdata,
    %  selecting this.PROC_FOLDER, converting to NIfTI, constructing json, and
    %  providing QC functionality.  Start with batch().
    %  
    %  Created 10-Feb-2022 19:47:37 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1837725 (R2021b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function [j,c] = parcluster()
            %% #PARCLUSTER
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
            %      subjectsDir (folder):  default ${ADNI_HOME}/loni_FDG/loni
            %      bidsDir (folder):  default ${ADNI_HOME}/bids/rawdata
            
            disp('mladni.AdniBidsFdg.batch starting')
            t0 = tic;

            subjectsDir = fullfile(getenv('ADNI_HOME'), 'loni_FDG', 'loni', '');
            bidsDir = fullfile(getenv('ADNI_HOME'), 'bids', 'rawdata', '');
            
            ip = inputParser;
            addParameter(ip, 'len', [], @isnumeric)
            addParameter(ip, 'subjectsDir', subjectsDir, @isfolder)
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:});
            ipr = ip.Results;      
            
            pwd0 = pushd(ipr.subjectsDir);
            globbed = glob(fullfile(ipr.subjectsDir, '*_S_*'))';
            if ~isempty(ipr.len)
                globbed = globbed(1:ipr.len);
            end
            fprintf('size(mladni.AdniBidsFdg.batch.globbed):  %s\n', mat2str(size(globbed)))
            
%            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)            
%                mladni.CHPC3.setenvs();
                procs = mladni.AdniBidsFdg.PROC_FOLDERS;
                for ip = 1:length(procs)
                    try
                        im = glob(fullfile(globbed{idx}, procs{ip}, '*', 'I*'))'; 
                        this = mladni.AdniBidsFdg();
                        for i = 1:length(im)
                            fn = this.build_nifti(im{i}, ipr.bidsDir); %#ok<PFBNS>
                            this.update_json(im{i}, ipr.bidsDir);
                            this.afni_3dresample(fn);
                        end
                    catch ME
                        handwarning(ME)
                    end
                end
            end
            popd(pwd0);
            
            t = toc(t0);
            disp('mladni.AdniBidsT1w.batch completed')  
        end
        function [tbl,csv_fn] = create_table_fdg_t1w(varargin)
            % Args:
            %      bidsDir (folder): default ${ADNI_HOME}/bids/rawdata
            % Returns:
            %      tbl: table for Matlab pairing fdg with best t1w for warping to atlases.
            %      csv_fn: exported for use in other platforms.

            ip = inputParser;
            bidsDir = fullfile(getenv('ADNI_HOME'), 'bids', 'rawdata', '');
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
        PROC_FOLDERS = {'Co-registered_Dynamic', 'Coreg,_Avg,_Std_Img_and_Vox_Siz,_Uniform_Resolution'}
    end

    properties (Dependent)
        table_pet_c3
        table_pet_qc
    end

    methods

        %% GET

        function g = get.table_pet_c3(this)
            if isempty(this.table_pet_c3_)
                this.table_pet_c3_ = this.adni_demographics.table_pet_c3();
            end
            g = this.table_pet_c3_;
        end
        function g = get.table_pet_qc(this)
            if isempty(this.table_pet_qc_)
                this.table_pet_qc_ = this.adni_demographics.table_pet_qc();
            end
            g = this.table_pet_qc_;
        end

        %%

        function fn = build_nifti(~, image_folder, rawdata_folder)
            %% calls dcm2niix
            %  Args:
            %      image_folder (folder): e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0/I12345/
            %      rawdata_folder (folder): e.g., /path/to/bids/rawdata
            %  Returns
            %      fn: e.g., /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-yyyyMMddHHmmss_trc-fdg_proc-CAS_pet.nii.gz

            assert(strcmp(image_folder(1), filesep)) % abs path

            % sub
            dt_folder = fileparts(fileparts(image_folder)); % e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0
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
            ensuredir(dest_folder);
            base = sprintf('%s_%s_trc-FDG_%s_pet', sub, ses14, proc);            
            if ~isfile(fullfile(dest_folder, [base '.nii.gz']))
                mlpipeline.Bids.dcm2niix(image_folder, 'f', base, 'o', dest_folder)
            end
            fn = fullfile(dest_folder, strcat(base, '.nii.gz'));
        end
        function p = quality(this, rid, dt)
            t = this.table_pet_qc(this.table_pet_qc.RID == rid & this.table_pet_qc.EXAMDATE == dt, :);
            if ~isempty(t)
                p = t.PASS;
                return
            end
            t = this.table_pet_c3(this.table_pet_c3.RID == rid & this.table_pet_c3.SCANDATE == dt, :);
            if ~isempty(t)
                p = t.SCANQLTY;
                return
            end
            p = nan;
        end
        function p = phase(this, rid, dt)
            t = this.table_pet_qc(this.table_pet_qc.RID == rid & this.table_pet_qc.EXAMDATE == dt, :);
            if ~isempty(t)
                p = t.Phase;
                return
            end
            t = this.table_pet_c3(this.table_pet_c3.RID == rid & this.table_pet_c3.SCANDATE == dt, :);
            if ~isempty(t)
                p = t.Phase;
                return
            end
            p = 'UNKNOWN';
        end
        function j = update_json(~, image_folder, rawdata_folder)
            %% calls jsonrecode
            %  Args:
            %      image_folder (folder): e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0/I12345/
            %      rawdata_folder (folder): e.g., /path/to/bids/rawdata
            %  Returns
            %      j: e.g., /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-yyyyMMddHHmmss_trc-fdg_proc-CAS_pet.nii.gz

            assert(strcmp(image_folder(1), filesep)) % abs path

            % sub
            dt_folder = fileparts(fileparts(image_folder)); % e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0
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
            s.ADNI_INFO.ScanDate = datestr(dt, 'yyyy-mm-dd'); %#ok<DATST> 
            s.ADNI_INFO.OriginalPath = image_folder;
            s.ADNI_INFO.OriginalFilename = basename(dcms{1});
            jsonrecode(j, s, 'filenameNew', j);
        end

        function this = AdniBidsFdg(varargin)
            this = this@mladni.AdniBidsT1w(varargin{:});
        end
    end
    
    %% PROTECTED

    properties (Access = protected)
        table_pet_qc_
        table_pet_c3_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
