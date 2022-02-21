classdef AdniBidsT1w < mladni.AdniBids
    %% line1
    %  line2
    %  
    %  Created 10-Feb-2022 19:47:50 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1837725 (R2021b) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function create(varargin)
            %% Creates, e.g., /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-20090713133153_acq-noaccel_proc-orig_T1s.json
            %                 /path/to/bids/rawdata/sub-128S0230/ses-20090713/pet/sub-128S0230_ses-20090713133153_trc-FDG_proc-CAS.nii.gz
            %  Also calls AdniBidsT1w.create_json_for_processed().
            %  Params:
            %      subjectsDir (folder):  e.g., $SINGULARITY_HOME/ADNI/loni_FDG/loni; default pwd.
            %      bidsDir (folder):  default $SINGULARITY_HOME?ADNI/bids/rawdata

            ip = inputParser;
            ip.KeepUnmatched = true;
            subjectsDir = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'loni_T1', 'loni', '');
            bidsDir = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'rawdata', '');
            addParameter(ip, 'subjectsDir', subjectsDir, @isfolder)
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            pwd0 = pushd(ipr.subjectsDir);
            this = mladni.AdniBidsT1w();
            s = globFoldersT(fullfile(ipr.subjectsDir, '*_S_*'));
            for si = 1:length(s)
                im = globFoldersT(fullfile(s{si}, '*/*-*-*_*_*_*/I*'));
                for i = 1:length(im)
                    if ~isempty(glob(fullfile(im{i}, '*.dcm')))
                        this.build_nifti(im{i}, ipr.bidsDir);
                    end
                    if ~isempty(glob(fullfile(im{i}, '*.nii')))
                        this.copy_nifti(im{i}, ipr.bidsDir);
                    end
                end
            end
            popd(pwd0);

            mladni.AdniBidsT1w.create_json_for_processed(varargin{:});
        end
        function create_json_for_processed(varargin)
            %  Params:
            %      bidsDir (folder):  default $SINGULARITY_HOME?ADNI/bids/rawdata

            import mladni.AdniBids.nii2json

            ip = inputParser;
            ip.KeepUnmatched = true;
            bidsDir = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'rawdata', '');
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            pwd0 = pushd(ipr.bidsDir);            
            anat = globFolders(fullfile(ipr.bidsDir, 'sub-*S*/ses-*/anat'));
            for ai = 1:length(anat)
                pwd1 = pushd(anat{ai});
                niis = globT('sub-*_ses-*_acq-*_proc-*.nii.gz');
                niis = niis(~contains(niis, '_proc-orig_'));
                for ni = 1:length(niis)
                    if ~isfile(nii2json(niis{ni}))
                        try
                            re = regexp(niis{ni}, 'sub-\d{3}S\d{4}_ses-\d{14}_acq-\w+_proc-(?<proc>\S+)_(T1w|IR).nii.gz', 'names');
                            nii_orig = strrep(niis{ni}, re.proc, 'orig');
                            copyfile(nii2json(nii_orig), nii2json(niis{ni}));
                        catch ME
                            handwarning(ME)
                        end
                    end
                end
                popd(pwd1);
            end
            popd(pwd0);
        end
    end

    methods
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

            % acq
            [~,seq_str] = fileparts(fileparts(dt_folder));
            if contains(seq_str, 'acc', 'IgnoreCase', true)
                acq = 'acq-accel';
            else
                acq = 'acq-noaccel';
            end

            % dcm2niix
            dest_folder = fullfile(rawdata_folder, sub, ses8, 'anat', '');
            if ~isfolder(dest_folder)
                mkdir(dest_folder);
            end
            cont = this.guess_contrast(seq_str);
            base = sprintf('%s_%s_%s_proc-orig_%s', sub, ses14, acq, cont);
            if ~isfile(fullfile(dest_folder, strcat(base, '.nii.gz')))
                s = []; r = '';
                try
                    [s,r] = mlpipeline.Bids.dcm2niix(image_folder, 'f', base, 'o', dest_folder);
                    this.update_json(image_folder, rawdata_folder);
                catch ME
                    if 0 ~= s
                        warning('mladni:RuntimeWarning', r)
                        handwarning(ME)
                    end
                end
            end                
        end
        function fn = copy_nifti(this, image_folder, rawdata_folder)
            %% copies nii from ADNI
            %  Args:
            %      image_folder (folder): e.g., /path/to/002_S_0413/{MP-RAGE,}/yyyy-MM-dd_HH_mm_ss.0/I12345/
            %      rawdata_folder (folder): e.g., /path/to/bids/rawdata
            %  Returns
            %      fn: e.g.,
            %      /path/to/bids/rawdata/sub-002S0413/ses-yyyyMMdd/sub-002S0413_ses-yyyy-MM-dd_HH_mm_ss_acq-accel_proc-gradbias_T1w.{json,nii.gz}

            assert(strcmp(image_folder(1), filesep)) % abs path

            % sub
            dt_folder = fileparts(image_folder);
            [~,sub_str] = fileparts(fileparts(fileparts(dt_folder)));
            sub = ['sub-' strrep(sub_str, '_', '')];

            % ses
            [~,dt_str] = fileparts(dt_folder);
            dt = datetime(strtok(dt_str, '.'), 'InputFormat', 'yyyy-MM-dd_HH_mm_ss');
            ses = sprintf('ses-%i%02i%02i%02i%02i%02i', ...
                year(dt), month(dt), day(dt), hour(dt), minute(dt), second(dt));
            ses8 = sprintf('ses-%i%02i%02i', ...
                year(dt), month(dt), day(dt));
            ses14 = sprintf('ses-%i%02i%02i%02i%02i%02i', ...
                year(dt), month(dt), day(dt), hour(dt), minute(dt), second(dt));

            % acq
            [~,seq_str] = fileparts(fileparts(dt_folder));
            if contains(seq_str, 'acc', 'IgnoreCase', true)
                acq = 'acq-accel';
            else
                acq = 'acq-noaccel';
            end

            % nii from ADNI, proc, cont
            dest_folder = fullfile(rawdata_folder, sub, ses8, 'anat', '');
            if ~isfolder(dest_folder)
                mkdir(dest_folder);
            end
            nii = globT(fullfile(image_folder, '*.nii'));
            fn = '';
            for i = 1:length(nii)
                proc = this.guess_proc(nii{i}, sub_str);
                cont = this.guess_contrast(nii{i});
                fn = fullfile(dest_folder, sprintf('%s_%s_%s_%s_%s.nii', sub, ses, acq, proc, cont));
                if ~isfile(strcat(fn, '.gz'))
                    copyfile(nii{1}, fn, 'f');
                    gzip(fn);
                    delete(fn);
                    fn = strcat(fn, '.gz');
                end
            end
        end
        function cont = guess_contrast(~, nii)
            if contains(nii, 'IR', 'IgnoreCase', true) && contains(nii, 'SPGR', 'IgnoreCase', true)
                cont = 'IR';
            else
                cont = 'T1w';
            end
        end
        function suff = guess_proc(~, nii, sub_str)
            assert(istext(nii))
            assert(istext(sub_str))
            suff = '';
            if contains(nii, 'mt1', 'IgnoreCase', true)
                suff = [suff '-mt1'];
            end
            if contains(nii, 'gradwarp', 'IgnoreCase', true)
                suff = [suff '-gradwarp'];
            end
            if contains(nii, 'n3', 'IgnoreCase', true)
                suff = [suff '-n3'];
            end
            if contains(nii, 'b1', 'IgnoreCase', true) && contains(nii, 'correction', 'IgnoreCase', true)
                suff = [suff '-b1corr'];
            end
            if contains(nii, 'scaled', 'IgnoreCase', true)
                suff = [suff '-scaled'];
            end
            if contains(nii, 'mask', 'IgnoreCase', true)
                suff = [suff '-mask'];
            end
            if contains(nii, 'regist', 'IgnoreCase', true)
                suff = [suff '-reg'];
            end

            % KLUDGE
            if isempty(suff)
                [~,fp_] = myfileparts(nii);
                suff = extractAfter(fp_, sub_str);
            end
            if isempty(suff)
                suff = 'unknown';
            end
            suff = strrep(suff, '_', '-');
            suff = ['proc' suff];
        end         
        function j = update_json(this, image_folder, rawdata_folder)
            %% calls jsonrecode
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

            % acq
            [~,seq_str] = fileparts(fileparts(dt_folder));
            if contains(seq_str, 'acc', 'IgnoreCase', true)
                acq = 'acq-accel';
            else
                acq = 'acq-noaccel';
            end

            % jsonrecode
            dest_folder = fullfile(rawdata_folder, sub, ses8, 'anat', '');
            cont = this.guess_contrast(seq_str);
            base = sprintf('%s_%s_%s_proc-orig_%s', sub, ses14, acq, cont);
            j = fullfile(dest_folder, strcat(base, '.json'));
            dcms = glob(fullfile(image_folder, '*.dcm'));
            try
                s.ADNI_INFO.ScanDate = datestr(dt, 'yyyy-mm-dd');
                s.ADNI_INFO.OriginalPath = image_folder;
                s.ADNI_INFO.OriginalFilename = basename(dcms{1});
                jsonrecode(j, s, 'filenameNew', j);                          
            catch
            end
        end

        function this = AdniBidsT1w(varargin)
            this = this@mladni.AdniBids(varargin{:})
            
            this.ad_ = mladni.AdniDemographics();
        end
    end
    
    %% PROTECTED

    properties (Access = protected)
        ad_
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
