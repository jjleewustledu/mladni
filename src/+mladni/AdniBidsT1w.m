classdef AdniBidsT1w < handle & mladni.AdniBids
    %% AdniBidsT1w migrates imaging from loni_T1 folders to bids/rawdata,
    %  converting to NIfTI, constructing json, and providing QC functionality.  Start with batch().
    %  
    %  Created 10-Feb-2022 19:47:50 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1837725 (R2021b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Dependent)
        table_mri_imageqc
        table_mri_quality
        table_mri_quality_adni3
    end

    methods

        %% GET

        function g = get.table_mri_imageqc(this)
            if isempty(this.table_mri_imageqc_)
                this.table_mri_imageqc_ = this.adni_merge.table_mri_imageqc();
            end
            g = this.table_mri_imageqc_;
        end
        function g = get.table_mri_quality(this)
            if isempty(this.table_mri_quality_)
                this.table_mri_quality_ = this.adni_merge.table_mri_quality();
            end
            g = this.table_mri_quality_;
        end
        function g = get.table_mri_quality_adni3(this)
            if isempty(this.table_mri_quality_adni3_)
                this.table_mri_quality_adni3_ = this.adni_merge.table_mri_quality_adni3();
            end
            g = this.table_mri_quality_adni3_;
        end

        %%

        function fn = build_nifti(this, image_folder, rawdata_folder)
            %% calls dcm2niix
            %  Args:
            %      image_folder (folder): e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0/I12345/
            %      rawdata_folder (folder): e.g., /path/to/bids/rawdata
            %  Returns
            %      fn: e.g.,
            %      /path/to/bids/rawdata/sub-002S0413/ses-yyyyMMdd/sub-002S0413_ses-yyyy-MM-dd_HH_mm_ss_acq-accel_proc-orig_T1w.{json,nii.gz}

            assert(strcmp(image_folder(1), filesep)) % abs path

            % sub
            dt_folder = fileparts(image_folder);
            [~,sub_str] = fileparts(fileparts(fileparts(dt_folder)));
            sub = ['sub-' strrep(sub_str, '_', '')];

            % ses
            [~,dt_str] = fileparts(dt_folder);
            dt = datetime(strtok(dt_str, '.'), 'InputFormat', 'yyyy-MM-dd_HH_mm_ss');
            ses8 = sprintf('ses-%i%02i%02i', ...
                year(dt), month(dt), day(dt));
            ses14 = sprintf('ses-%i%02i%02i%02i%02i%02i', ...
                year(dt), month(dt), day(dt), hour(dt), minute(dt), second(dt));

            % dcm2niix
            dest_folder = fullfile(rawdata_folder, sub, ses8, 'anat', '');
            ensuredir(dest_folder);
            [~,seq_str] = fileparts(fileparts(dt_folder));
            g = glob(fullfile(image_folder, '*.dcm'));
            proc_str = mybasename(g{1});
            acq = this.guess_acq(strcat(seq_str, proc_str));
            proc = this.guess_proc(proc_str, sub_str);
            if contains(proc, '-cal')
                fn = '';
                return
            end
            cont = this.guess_contrast(seq_str);
            if strcmp(cont, 'local')
                fn = '';
                return
            end
            base = sprintf('%s_%s_%s_%s-dcm_%s', sub, ses14, acq, proc, cont);
            baserpi = sprintf('%s_%s_%s_%s-dcm_orient-rpi_%s', sub, ses14, acq, proc, cont);
            if ~isfile(fullfile(dest_folder, strcat(base, '.nii.gz'))) && ...
               ~isfile(fullfile(dest_folder, strcat(baserpi, '.nii.gz')))
                s = []; r = '';
                try
                    [s,r] = mlpipeline.Bids.dcm2niix(image_folder, 'f', base, 'o', dest_folder, 'version', 20180627);
                    fn = this.update_json(image_folder, rawdata_folder, base);
                    fn = strrep(fn, '.json', '.nii.gz');
                    assert(0 == s);
                catch ME
                    if 0 ~= s
                        [s,r] = mlpipeline.Bids.dcm2niix(image_folder, 'f', base, 'o', dest_folder);   
                    end
                    if 0 ~= s
                        warning('mladni:RuntimeWarning', r)
                        handwarning(ME)
                    end
                end
            end                
        end
        function fngz = copy_nifti(this, image_folder, rawdata_folder, opts)
            %% copies nii from ADNI
            %  Args:
            %      this mladni.AdniBids
            %      image_folder {mustBeFolder} : source
            %      rawdata_folder {mustBeFolder} : destination
            %      opts.noclobber logical = false : in particular, don't clobber .json
            %  Returns
            %      fngz: e.g.,
            %      /path/to/bids/rawdata/sub-002S0413/ses-yyyyMMdd/sub-002S0413_ses-yyyy-MM-dd_HH_mm_ss_acq-accel_proc-gradbias_T1w.{json,nii.gz}

            arguments
                this mladni.AdniBids
                image_folder {mustBeFolder}
                rawdata_folder {mustBeFolder}
                opts.noclobber logical = false
            end
            assert(strcmp(image_folder(1), filesep)) % abs path

            % sub
            dt_folder = fileparts(image_folder);
            [~,sub_str] = fileparts(fileparts(fileparts(dt_folder)));
            sub = ['sub-' strrep(sub_str, '_', '')];

            % ses
            [~,dt_str] = fileparts(dt_folder);
            dt = datetime(strtok(dt_str, '.'), 'InputFormat', 'yyyy-MM-dd_HH_mm_ss');
            ses14 = sprintf('ses-%i%02i%02i%02i%02i%02i', ...
                year(dt), month(dt), day(dt), hour(dt), minute(dt), second(dt));
            ses8 = sprintf('ses-%i%02i%02i', ...
                year(dt), month(dt), day(dt));
            %ses14 = sprintf('ses-%i%02i%02i%02i%02i%02i', ...
            %    year(dt), month(dt), day(dt), hour(dt), minute(dt), second(dt));

            % nii from ADNI, proc, cont
            dest_folder = fullfile(rawdata_folder, sub, ses8, 'anat', '');
            ensuredir(dest_folder);
            [~,seq_str] = fileparts(fileparts(dt_folder));
            nii = globT(fullfile(image_folder, '*.nii'));
            fn = '';
            for i = 1:length(nii)
                acq = this.guess_acq(strcat(seq_str, nii{i}));
                proc = this.guess_proc(nii{i}, sub_str);
                if contains(proc, '-mask') || contains(proc, '-masked')
                    continue
                end
                if contains(proc, '-epi')
                    continue
                end
                if contains(proc, '-cal')
                    continue
                end
                cont = this.guess_contrast(nii{i});
                if strcmp(cont, 'local')
                    continue
                end
                fn = fullfile(dest_folder, sprintf('%s_%s_%s_%s_%s.nii', sub, ses14, acq, proc, cont));
                fngz = strcat(fn, '.gz');
                fnrpi = fullfile(dest_folder, sprintf('%s_%s_%s_%s_orient-rpi_%s.nii', sub, ses14, acq, proc, cont));
                fnrpigz = strcat(fnrpi, '.gz');

                if ~isfile(fngz) && ~isfile(fnrpigz)
                    copyfile(nii{1}, fn, 'f');
                    fngz = ensuregz(fn);
                end

                % de novo json, recording ADNI_INFO for processed NIfTI from LONI
                fnj = strrep(fngz, '.nii.gz', '.json');
                fnrpij = strrep(fnrpigz, '.nii.gz', '.json');
                if ~opts.noclobber || (~isfile(fnj) && ~isfile(fnrpij))
                   s.ADNI_INFO.ScanDate = char(datetime(dt, Format='yyyy-MM-dd'));
                   s.ADNI_INFO.OriginalPath = myfileparts(nii{1});
                   s.ADNI_INFO.OriginalFilename = mybasename(strcat(nii{1}, '.nii'));
                   fid = fopen(fnj, 'w');
                   fprintf(fid, jsonencode(s, 'PrettyPrint', true));
                   fclose(fid);
                end
            end
        end
        function acq = guess_acq(~, str)
            acq = 'acq-unk';
            if contains(str, '_Acc', 'IgnoreCase', true)
                acq = 'acq-accel';
            end
            if contains(str, '_SENSE', 'IgnoreCase', true)
                acq = 'acq-sense';
            end
            if contains(str, '_GRAPPA', 'IgnoreCase', true)
                acq = 'acq-grappa';
            end
        end
        function cont = guess_contrast(~, nii)
            if contains(nii, 'Localiz', 'IgnoreCase', true)
                cont = 'local';
                return
            end
            if contains(nii, 'IR', 'IgnoreCase', true) && contains(nii, 'SPGR', 'IgnoreCase', true)
                cont = 'irspgr_T1w';
                return
            end
            if contains(nii, 'MT1', 'IgnoreCase', true)
                cont = 'mt_T1w';
                return
            end
            if contains(nii, 'T2', 'IgnoreCase', true)
                cont = 'T2w';
                return
            end
            cont = 'T1w';            
        end
        function suff = guess_proc(~, nii, sub_str)
            assert(istext(nii))
            assert(istext(sub_str))
            suff = '';
            if contains(nii, 'b1', 'IgnoreCase', true) && contains(nii, 'correction', 'IgnoreCase', true)
                suff = [suff '-b1corr'];
            end
            %if contains(nii, '_Br_', 'IgnoreCase', true)
            %    suff = [suff '-br'];
            %end
            if contains(nii, 'Calibration', 'IgnoreCase', true)
                suff = [suff '-cal'];
            end
            if contains(nii, '_EPI_', 'IgnoreCase', true)
                suff = [suff '-epi'];
            end
            if contains(nii, 'gradwarp', 'IgnoreCase', true)
                suff = [suff '-gradwarp'];
            end
            %if contains(nii, '_ICS_', 'IgnoreCase', true)
            %    suff = [suff '-ics'];
            %end
            if contains(nii, 'mask', 'IgnoreCase', true)
                suff = [suff '-mask'];
            end
            if contains(nii, 'mt1', 'IgnoreCase', true)
                suff = [suff '-mt1'];
            end
            if contains(nii, 'n3', 'IgnoreCase', true)
                suff = [suff '-n3'];
            end
            if contains(nii, 'n4', 'IgnoreCase', true)
                suff = [suff '-n4'];
            end
            if contains(nii, 'scaled', 'IgnoreCase', true)
                suff = [suff '-scaled'];
            end
            if contains(nii, 'spatially_normalized,_masked', 'IgnoreCase', true)
                suff = [suff '-masked'];
            end
            if contains(nii, '_raw_', 'IgnoreCase', true)
                suff = [suff '-raw'];
            end
            if contains(nii, 'regist', 'IgnoreCase', true)
                suff = [suff '-reg'];
            end
            if contains(nii, 'REPE', 'IgnoreCase', true)
                suff = [suff '-rep'];
            end

            % KLUDGE
            if isempty(suff)
                [~,fp_] = myfileparts(nii);
                suff = extractAfter(fp_, sub_str);
                suff = strrep(suff, '_', '-');
            end
            if isempty(suff)
                suff = '-unknown';
            end
            suff = ['proc' suff];
        end         
        function fn = select_t1w(this, t1w_arr, fdg_fn)
            %  Args:
            %      t1w_arr cell: t1w fqfilenames for a subject
            %      fdg_fn(file): fdg fqfilename for a subject
            %  Returns:
            %      fn: selected from t1w_arr

            arguments
                this mladni.AdniBidsT1w 
                t1w_arr cell
                fdg_fn {mustBeTextScalar}
            end
            fdg_dt = this.filestr2datetime(basename(fdg_fn));

            % prune
            [t1w_arr,t1w_def] = this.select_t1w_by_quality(t1w_arr); 
            if isempty(t1w_arr)
                t1w_arr = t1w_def; 
            end
            [t1w_arr,dur] = this.sort_t1w_by_deltadate(t1w_arr, fdg_dt);
            assert(~isempty(dur), stackstr(2))
            dur = hours(dur);

            % rank-order T1w by repeat scans, acq. accel., n3, not magn. transfer
            % e.g.:  acq_sense preferred over acq_unk
            rep = cell2mat(cellfun(@(x) contains(x, '-rep'), t1w_arr, UniformOutput=false));
            acq = cellfun(@(x) regexp(x, 'acq-[a-z]+', 'match'), t1w_arr, UniformOutput=false);
            acq = cellfun(@(x) length(x{1}), acq, UniformOutput=false);
            acq = ascol([acq{:}]);
            n3 = cell2mat(cellfun(@(x) contains(x, '-n3'), t1w_arr, UniformOutput=false));
            mt1 = cell2mat(cellfun(@(x) contains(x, '-mt1'), t1w_arr, UniformOutput=false));
            t_ = table(t1w_arr, dur, acq, n3, mt1, rep);

            % Failure on R2021a.
            % Warning: An error occurred when sorting rows of the variable 'rep', whose type is 'cell'.
            % (Type "warning off MATLAB:table:sortrows:SortrowsOnVarFailed" to suppress this warning.)            
            t_ = sortrows(t_, {'dur', 'acq', 'n3', 'mt1', 'rep'}, {'ascend' 'descend' 'descend' 'ascend' 'descend'});

            fn = t_.t1w_arr{1};
        end
        function [t1w_arr,t1w_def] = select_t1w_by_quality(this, t1w_arr)
            %% selects PASS == 1, and rejects tags for mask, epi, cal, local.
            %  Args:
            %      t1w_arr cell: t1w fqfilenames for subject(s)
            %  Returns:
            %      t1w_arr cell: select PASS
            %      t1w_def cell: select ~PASS

            arguments
                this mladni.AdniBidsT1w
                t1w_arr cell
            end
            if isempty(t1w_arr)
                return
            end

            t1w_arr = t1w_arr(~contains(t1w_arr, '-mask')); % for safety
            t1w_arr = t1w_arr(~contains(t1w_arr, '-epi')); % for safety
            t1w_arr = t1w_arr(~contains(t1w_arr, '-cal')); % for safety
            t1w_arr = t1w_arr(~contains(t1w_arr, 'local')); % for safety
            t1w_arr = t1w_arr(~contains(t1w_arr, 'Coreg')); % for safety
            t1w_arr = t1w_arr(~contains(t1w_arr, 'Co-reg')); % for safety
            t1w_def = [];

            return

            %% PASS, STUDY_OVERALLPASS
            rids = cell2mat(cellfun(@this.filestr2rid, t1w_arr, UniformOutput=false));
            rids = unique(rids);
            selected = ismember(this.table_mri_quality.RID, rids);
            q = this.table_mri_quality(selected, :);
            selected3 = ismember(this.table_mri_quality_adni3.RID, rids);
            q3 = this.table_mri_quality_adni3(selected3, :);
            if all(q.PASS == 1) && all(q3.STUDY_OVERALLPASS == 1)
                return
            end

            defects = false(size(t1w_arr));
            q = q(q.PASS ~= 1, :);
            for qi = 1:size(q,1)
                defects = defects | ...
                    (this.samedate(this.cell2datetimes(t1w_arr), q.EXAMDATE(qi)) & ...
                     strcmp(this.fqcell2uid(t1w_arr), q.LONIUID(qi)));
            end
            q3 = q3(q3.STUDY_OVERALLPASS ~= 1, :);
            for qi = 1:size(q3,1)
                defects = defects | ...
                    this.samedate(this.cell2datetimes(t1w_arr), q3.SERIES_DATE(qi));
            end
            assert(~isempty(t1w_arr), stackstr(2))

            t1w_arr = t1w_arr(~defects);
            t1w_def = t1w_arr(defects);
        end
        function [t1w_arr,dur] = sort_t1w_by_deltadate(this, t1w_arr, dt)
            %  Args:
            %      t1w_arr cell: t1w fqfilenames for subject(s)
            %      dt datetime = NaT: datetime of associated event, e.g., FDG scan
            %  Returns:
            %      t1w_arr: nearest to dt at 1st element
            %      dur: absolute durations of t1w_arr from dt

            arguments
                this mladni.AdniBidsT1w 
                t1w_arr cell
                dt datetime = NaT
            end
            if isempty(t1w_arr) || isnat(dt)
                dur = [];
                return
            end
            
            try
                t1w_dt = ascol(this.cell2datetimes(t1w_arr));
                dur = ascol(abs(t1w_dt - dt));
                t_ = table(t1w_arr, dur);
                [t_,idx] = sortrows(t_, 'dur', 'ascend');
                t1w_arr = t_.t1w_arr;
                dur = dur(idx);
            catch ME
                handwarning(ME)
            end
        end
        function t1w_arr = sort_t1w_by_proc(this, t1w_arr)
            %  Args:
            %      t1w_arr cell: t1w fqfilenames for subject(s)
            %  Returns:
            %      t1w_arr: rep, length(acq), n3, ~mt1 at 1st element

            arguments
                this mladni.AdniBidsT1w %#ok<INUSA> 
                t1w_arr cell
            end

            % rank-order T1w by repeat scans, acq. accel., n3, not magn. transfer
            rep = cellfun(@(x) contains(x, '-rep'), t1w_arr, UniformOutput=false);
            rep = ascol([rep{:}]);
            acq = cellfun(@(x) regexp(x, 'acq-[a-z]+', 'match'), t1w_arr, UniformOutput=false);
            acq = cellfun(@(x) length(x{1}), acq, UniformOutput=false);
            acq = ascol([acq{:}]);
            n3 = cellfun(@(x) contains(x, '-n3'), t1w_arr, UniformOutput=false);
            n3 = ascol([n3{:}]);
            mt1 = cellfun(@(x) contains(x, '-mt1'), t1w_arr, UniformOutput=false);
            mt1 = ascol([mt1{:}]);
            t_ = table(t1w_arr, rep, acq, n3, mt1);
            t_ = sortrows(t_, {'rep', 'acq', 'n3', 'mt1'}, {'descend' 'descend' 'descend' 'ascend'});

            t1w_arr = t_(:,1);
        end
        function j = update_json(~, image_folder, rawdata_folder, base)
            %% calls jsonrecode
            %  Args:
            %      image_folder (folder): e.g., /path/to/002_S_0413/Co*/yyyy-MM-dd_HH_mm_ss.0/I12345/
            %      rawdata_folder (folder): e.g., /path/to/bids/rawdata
            %      base (text): e.g., sub-002S0295_ses-20060418082030_acq-noaccel_proc-unknown-dcm_T1w
            %  Returns
            %      fn: e.g.,
            %      /path/to/bids/rawdata/sub-002S0413/ses-yyyyMMdd/sub-002S0413_ses-yyyy-MM-dd_HH_mm_ss_acq-accel_proc-orig_T1w.{json,nii.gz}

            assert(strcmp(image_folder(1), filesep)) % abs path

            % sub
            dt_folder = fileparts(image_folder);
            [~,sub_str] = fileparts(fileparts(fileparts(dt_folder)));
            sub = ['sub-' strrep(sub_str, '_', '')];

            % ses
            [~,dt_str] = fileparts(dt_folder);
            dt = datetime(strtok(dt_str, '.'), 'InputFormat', 'yyyy-MM-dd_HH_mm_ss');
            ses8 = sprintf('ses-%i%02i%02i', ...
                year(dt), month(dt), day(dt));

            % jsonrecode
            dest_folder = fullfile(rawdata_folder, sub, ses8, 'anat', '');
            j = fullfile(dest_folder, strcat(base, '.json'));
            dcms = glob(fullfile(image_folder, '*.dcm'));
            try
                s.ADNI_INFO.ScanDate = datestr(dt, 'yyyy-mm-dd'); %#ok<DATST> 
                s.ADNI_INFO.OriginalPath = image_folder;
                s.ADNI_INFO.OriginalFilename = basename(dcms{1});
                jsonrecode(j, s, 'filenameNew', j);                          
            catch
            end
        end

        function this = AdniBidsT1w(varargin)
            this = this@mladni.AdniBids(varargin{:});
        end
    end
    
    methods (Static)
        function propcluster()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '6000'; % in MB; deepmrseg requires 10 GB; else 6 GB
            c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = '';
            c.AdditionalProperties.WallTime = '72:00:00'; % 24 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function [j,c] = parcluster()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/          
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.AdniBidsT1w.batch, 1, {'len', []}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false); % {'len', []}, 'Pool', 31            
        end
        function [j,c] = parcluster2()
            %% #PARCLUSTER2
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/          
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.AdniBidsT1w.batch2, 1, {'len', []}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false); % {'len', []}, 'Pool', 31            
        end
        function t = batch(varargin)
            %% From subjectsDir, dcm -> nifti or nifti -> nifti, for ADNI preprocessed
            %  Params:
            %      len (numeric):  # subjects
            %      subjectsDir (folder):  default $ADNI_HOME/loni_T1/loni
            %      bidsDir (folder):  default ${ADNI_HOME}/bids/rawdata
            
            disp('mladni.AdniBidsT1w.batch starting')
            t0 = tic;

            mladni.CHPC3.setenvs();
            subjectsDir = fullfile(getenv('ADNI_HOME'), 'loni_T1', 'loni', '');
            bidsDir = fullfile(getenv('ADNI_HOME'), 'bids', 'rawdata', '');
            
            ip = inputParser;
            addParameter(ip, 'len', [], @isnumeric)
            addParameter(ip, 'subjectsDir', subjectsDir, @isfolder)
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:});
            ipr = ip.Results;      
            
            pwd0 = pushd(ipr.subjectsDir);
            %select = {'067_S_4782' '067_S_4918' '067_S_5159' '067_S_5160' '067_S_5205' ...
            %          '067_S_5212' '067_S_6474' '067_S_6525' '067_S_6529' '068_S_0127'};
            %globbed = cellfun(@(x) fullfile(ipr.subjectsDir, x), select, 'UniformOutput', false);            
            globbed = globFoldersT(fullfile(ipr.subjectsDir, '*_S_*'));
            if ~isempty(ipr.len)
                globbed = globbed(1:ipr.len);
            end
            fprintf('size(mladni.AdniBidsT1w.batch.globbed):  %s\n', mat2str(size(globbed)))
            
            %mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)
                mladni.CHPC3.setenvs();
                im = globFoldersT(fullfile(globbed{idx}, '*/*-*-*_*_*_*/I*'));                    
                this = mladni.AdniBidsT1w();
                for i = 1:length(im)
                    try
                        if contains(im{i}, 'localizer', 'IgnoreCase', true)
                            continue
                        end
                        if contains(im{i}, 'Coreg', 'IgnoreCase', true) || ...
                           contains(im{i}, 'Co-reg', 'IgnoreCase', true)
                            continue
                        end
                        if ~isempty(glob(fullfile(im{i}, '*.dcm')))
                            fn = this.build_nifti(im{i}, ipr.bidsDir); %#ok<PFBNS>
                            mladni.AdniBids.afni_3dresample(fn);
                        end
                        if ~isempty(glob(fullfile(im{i}, '*.nii')))
                            fn = this.copy_nifti(im{i}, ipr.bidsDir);
                            mladni.AdniBids.afni_3dresample(fn);
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
        function t = batch2(varargin)
            %% Requires completion of batch(). Generates json for LONI preprocessed files by searching for matching 
            %  original files.
            %  Params:
            %      len (numeric):  # sessions
            %      bidsDir (folder):  default ${ADNI_HOME}/bids/rawdata

            import mladni.AdniBids.nii2json    
            disp('Start mladni.AdniBidsT1w.batch2()')            
            t0 = tic;

            mladni.CHPC3.setenvs();
            bidsDir = fullfile(getenv('ADNI_HOME'), 'bids', 'rawdata_20221202', '');
            
            ip = inputParser;
            addParameter(ip, 'len', [], @isnumeric)
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            pwd0 = pushd(ipr.bidsDir);            
            globbed = globFolders(fullfile(ipr.bidsDir, 'sub-*S*/ses-*/anat'));
            if ~isempty(ipr.len)
                globbed = globbed(1:ipr.len);
            end
            fprintf('mladni.AdniBidsT1w.batch2.globbed:\n')
            disp(size(globbed))
            
            mladni.CHPC3.clean_tempdir();
            for idx = 1:length(globbed)
                mladni.CHPC3.setenvs();
                pwd1 = pushd(globbed{idx});
                niis = globT('sub-*_ses-*_acq-*_proc-*.nii.gz');
                for ni = 1:length(niis)
                    if ~isfile(nii2json(niis{ni}))
                        try
                            re = regexp(niis{ni}, '(?<sub>sub-\d{3}S\d{4})_(?<ses>ses-\d{14})_\S+.nii.gz', 'names');
                            globbed_json = glob(sprintf('%s_%s_*.json', re.sub, re.ses));
                            if isempty(globbed_json)
                                continue
                            end
                            mlbash(sprintf('cp -f %s %s', globbed_json{end}, nii2json(niis{ni})));
                        catch ME
                            handwarning(ME)
                        end
                    end
                end
                popd(pwd1);
            end
            popd(pwd0);
            
            t = toc(t0);
            disp('mladni.AdniBidsT1w.batch2() completed')  
        end
        function create(varargin)
            %% Creates, e.g., /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-20090713133153_acq-noaccel_proc-orig_T1s.json
            %                 /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-20090713133153_trc-FDG_proc-CAS.nii.gz
            %  Also calls AdniBidsT1w.create_json_for_processed().
            %  Params:
            %      subjectsDir (folder):  e.g., ${ADNI_HOME}/loni_FDG/loni; default pwd.
            %      bidsDir (folder):  default ${ADNI_HOME}/bids/rawdata

            ip = inputParser;
            ip.KeepUnmatched = true;
            subjectsDir = fullfile(getenv('ADNI_HOME'), 'loni_T1', 'loni', '');
            bidsDir = fullfile(getenv('ADNI_HOME'), 'bids', 'rawdata', '');
            addParameter(ip, 'subjectsDir', subjectsDir, @isfolder)
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            pwd0 = pushd(ipr.subjectsDir);
            s = globFoldersT(fullfile(ipr.subjectsDir, '*_S_*'));
            parfor si = 1:length(s)
                im = globFoldersT(fullfile(s{si}, '*/*-*-*_*_*_*/I*'));
                for i = 1:length(im)
                    if contains(im{i}, 'localizer', 'IgnoreCase', true)
                        continue
                    end
                    if contains(im{i}, 'Coreg', 'IgnoreCase', true) || ...
                       contains(im{i}, 'Co-reg', 'IgnoreCase', true)
                        continue
                    end
                    if ~isempty(glob(fullfile(im{i}, '*.dcm')))
                        this = mladni.AdniBidsT1w();
                        this.build_nifti(im{i}, ipr.bidsDir);
                    end
                    if ~isempty(glob(fullfile(im{i}, '*.nii')))
                        this = mladni.AdniBidsT1w();
                        this.copy_nifti(im{i}, ipr.bidsDir);
                    end
                end
            end
            popd(pwd0);

            mladni.AdniBidsT1w.afni_3dresample_bidsdir(varargin{:});
        end
        function create_json_for_processed(varargin)
            %% Requires completion of copy_nifti().
            %  Params:
            %      bidsDir (folder):  default ${ADNI_HOME}/bids/rawdata

            import mladni.AdniBids.nii2json

            ip = inputParser;
            ip.KeepUnmatched = true;
            bidsDir = fullfile(getenv('ADNI_HOME'), 'bids', 'rawdata', '');
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            pwd0 = pushd(ipr.bidsDir);            
            anat = globFolders(fullfile(ipr.bidsDir, 'sub-*S*/ses-*/anat'));
            for ai = 1:length(anat)
                pwd1 = pushd(anat{ai});
                niis = globT('sub-*_ses-*_acq-*_proc-*.nii.gz');
                niis = niis(~contains(niis, '_proc-dcm'));
                for ni = 1:length(niis)
                    if ~isfile(nii2json(niis{ni})) % missing json, likely b/c copied from NIfTI
                        try
                            re = regexp(niis{ni}, '(?<subses>sub-\d{3}S\d{4}_ses-\d{14})\S+.nii.gz', 'names');
                            json_orig = glob(sprintf(strcat(re.subses, '*.json')));
                            copyfile(json_orig{end}, nii2json(niis{ni}));
                        catch ME
                            handwarning(ME)
                        end
                    end
                end
                popd(pwd1);
            end
            popd(pwd0);
        end
        function afni_3dresample_bidsdir(varargin)
            %  Params:
            %      bidsDir (folder):  default ${ADNI_HOME}/bids/rawdata

            import mladni.AdniBids.nii2json

            ip = inputParser;
            ip.KeepUnmatched = true;
            bidsDir = fullfile(getenv('ADNI_HOME'), 'bids', 'rawdata', '');
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            pwd0 = pushd(ipr.bidsDir);            
            anat = globFolders(fullfile(ipr.bidsDir, 'sub-*S*/ses-*/anat'));
            for ai = 1:length(anat)
                pwd1 = pushd(anat{ai});
                niis = globT('sub-*_ses-*.nii.gz');
                for ni = 1:length(niis)
                    if ~contains(niis{ni}, '_orient-rpi')
                        try
                            mladni.AdniBids.afni_3dresample(niis{ni});
                            deleteExisting(niis{ni});
                        catch ME
                            handwarning(ME)
                        end
                    end
                end

                % clean leftover json without nii.gz
                js = globT('sub-*_ses-*.json');
                for ji = 1:length(js)
                    its_nii = strrep(js{ji}, '.json', '.nii.gz');
                    if ~isfile(its_nii)
                        deleteExisting(js{ji})
                    end
                end

                popd(pwd1);
            end
            popd(pwd0);
        end
    end

    %% PROTECTED

    properties (Access = protected)
        table_mri_quality_
        table_mri_quality_adni3_
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
