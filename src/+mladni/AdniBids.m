classdef AdniBids < handle
    %% provides supporting methods for class hierarchy
    %  
    %  Created 30-Dec-2021 16:40:32 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1809720 (R2021b) Update 1 for MACI64.  Copyright 2021 John J. Lee.
    
    properties (Dependent)
        adni_demographics
        adni_merge
    end
   
    methods 

        %% GET

        function g = get.adni_demographics(this)
            if isempty(this.adni_demo_)
                this.adni_demo_ = mladni.AdniDemographics();
            end
            g = this.adni_demo_;
        end
        function g = get.adni_merge(this)
            g = this.adni_demographics.adni_merge;
        end
    end

    methods (Static)
        % implementing AFNI 3dresample
        function [j,c] = parcluster_3dresample()
            %% #parcluster_3dresample
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/          
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.AdniBids.batch_3dresample, 1, {'len', []}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false); % {'len', []}, 'Pool', 31            
        end
        function t = batch_3dresample(varargin)
            %% Requires well-formed bidDir from batch() and batch2().  Applies 3dresample, updating bids filenames and
            %  json file contents.
            %  Params:
            %      len (numeric):  # sessions
            %      bidsDir (folder):  default ${ADNI_HOME}/bids/rawdata

            import mladni.AdniBids.nii2json
            disp('Start mladni.AdniBids.batch_3dresample()')            
            t0 = tic;
            
            bidsDir = fullfile(getenv('ADNI_HOME'), 'bids', 'rawdata', '');

            ip = inputParser;
            addParameter(ip, 'len', [], @isnumeric)
            addParameter(ip, 'bidsDir', bidsDir, @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;

            pwd0 = pushd(ipr.bidsDir);            
            globbed = globFolders(fullfile(ipr.bidsDir, 'sub-*S*/ses-*/*'));
            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)
                mladni.CHPC3.setenvs();
                pwd1 = pushd(globbed{idx});
                niis = globT('sub-*_ses-*.nii.gz');
                for ni = 1:length(niis)
                    if ~contains(niis{ni}, '_orient-rpi')
                        try
                            mladni.AdniBids.afni_3dresample(niis{ni});
                        catch ME
                            handwarning(ME)
                        end
                    end
                end
                popd(pwd1);
            end
            popd(pwd0);
            
            t = toc(t0);
            disp('mladni.AdniBids.batch_3dresample completed')  
        end        
        function fn1 = afni_3dresample(fn)
            fn1 = mlpipeline.Bids.afni_3dresample(fn);
        end

        % support methods
        function dt = cell2datetimes(c)
            if ischar(c) || isstring(c)
                c = {c};
            end
            dt = NaT(size(c));
            dt = ensureTimeZone(dt);
            for fi = 1:length(c)
                dt(fi) = mladni.AdniBids.filestr2datetime(c{fi}); 
            end
        end
        function dt = filestr2datetime(str)
            %% file or folder names containing 'ses-yyyyMMdd+' -> 'yyyyMMdd+'
            %  Args:
            %      str {mustBeTextScalar}

            dt = mladni.AdniBids.filestr2tagged('ses', str);
            switch length(dt)
                case 16:21
                    dt = datetime(dt, InputFormat='yyyyMMddHHmmss.S');
                case 14
                    dt = datetime(dt, InputFormat='yyyyMMddHHmmss');
                case 8
                    dt = datetime(dt, InputFormat='yyyyMMdd');
                otherwise
                    error('mladni:ValueError', stackstr(2))
            end
            dt = ensureTimeZone(dt);
        end
        function a = filestr2acq(str)
            %% file or folder names containing '_acq-grappa_' -> 'grappa'
            %  Args:
            %      str {mustBeTextScalar}

            a = mladni.AdniBids.filestr2tagged('acq', str);
        end
        function o = filestr2orient(str)
            %% file or folder names containing '_orient-rpi_' -> 'rpi'
            %  Args:
            %      str {mustBeTextScalar}

            o = mladni.AdniBids.filestr2tagged('orient', str);
        end
        function p = filestr2proc(str)
            %% file or folder names containing '_proc-aproc-anotherproc_' -> 'aproc-anotherproc'
            %  Args:
            %      str {mustBeTextScalar}

            p = mladni.AdniBids.filestr2tagged('proc', str);
        end
        function rid = filestr2rid(str)
            import mladni.AdniBids.subfold2rid
            import mladni.AdniBids.filestr2subfold
            rid = subfold2rid( ...
                filestr2subfold(str));
        end
        function s = filestr2ses(str)
            %% file or folder names containing '_ses-202201010000_' -> '202201010000'
            %  Args:
            %      str {mustBeTextScalar}

            s = mladni.AdniBids.filestr2tagged('ses', str);
        end
        function s = filestr2sub(str)
            %% file or folder names containing '_sub-012S3456_' -> '012S3456'
            %  Args:
            %      str {mustBeTextScalar}

            s = mladni.AdniBids.filestr2tagged('sub', str);
        end
        function fold = filestr2subfold(str)
            sub_item = mladni.AdniBids.filestr2sub(str);
            fold = strcat('sub-', sub_item);
        end
        function p = filestr2tagged(tag, str)
            %% file or folder names containing '_tag-anitem-anotheritem_' -> 'anitem-anotheritem'
            %  Args:
            %      tag {mustBeTextScalar}
            %      str {mustBeTextScalar}

            arguments
                tag {mustBeTextScalar}
                str {mustBeTextScalar}
            end
            ss = strsplit(str, filesep);
            ss = ss(contains(ss, strcat(tag, "-")));
            str = ss{end};

            re = regexp(str, sprintf("\\S*%s-(?<tagged>[a-zA-Z0-9\\-]+)_\\S+", tag), "names");
            p = re.tagged;
        end
        function uid = fqcell2uid(fqniis)
            if ischar(fqniis) || isstring(fqniis)
                fqniis = {fqniis};
            end
            uid = cell(size(fqniis));
            for fi = 1:length(fqniis)
                try
                    [pth,fp] = myfileparts(fqniis{fi});
                    j = fullfile(pth, strcat(fp, '.json'));
                    s = jsondecode(fileread(j));
                    re = regexp(s.ADNI_INFO.OriginalFilename, '\w+_(?<uid>S\d+)_I\d+.dcm', 'names');
                    uid{fi} = re.uid;
                catch ME
                    handwarning(ME)
                    uid{fi} = 'unknown';
                end
            end
        end
        function move_niis(varargin)
            %  Params:
            %      niis (cell): e.g., {'/path/to/loni/128_S_0230/sub_128S0230_ses-yyyyMMddHHmmss_trc-FDG_proc-CAS_pet.nii.gz' ...}
            %      dest (folder): e.g., '/path/to/adni/rawdata/sub-128S0230/ses-20120823/pet'
            %      table_sub (table): subject selected from this.table_fdg1_.

            import mladni.AdniBids.cell2datetimes
            import mladni.AdniBids.viscode2months

            ip = inputParser;
            addRequired(ip, 'niis', @iscell)
            addParameter(ip, 'dest', '', @isfolder)
            addParameter(ip, 'table_sub', [], @istable)
            addParameter(ip, 'viscode2', '', @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;

            if isempty(ipr.viscode2)
                    niis = ipr.niis;
                    for i = 1:length(niis)
                        try
                            movefile(niis{i}, ipr.dest, 'f')        
                        catch ME
                            handwarning(ME)
                            disp(niis{i})
                        end
                        try
                            [pth,fp] = myfileparts(niis{i});
                            j = fullfile(pth, strcat(fp, '.json'));
                            movefile(j, ipr.dest, 'f')
                        catch ME
                            handwarning(ME)
                            disp(j)
                        end
                    end
                return
            end

            try
                dt_bl = ipr.table_sub.AcqDate(contains(ipr.table_sub.VISCODE2, 'bl'));
                dts = cell2datetimes(ipr.niis);
                mm = calmonths(between(dt_bl, dts));
                mm(mm < 0) = 0;
                m_dest = viscode2months(ipr.viscode2);
                for m_delta = 1:40
                    niis = ipr.niis(abs(mm - m_dest) <= m_delta);
                    if isempty(niis)
                        continue
                    end
                end
                for i = 1:length(niis)
                    [~,fp,ext] = myfileparts(niis{i});
                    dest_fqfn = fullfile(ipr.dest, strcat(fp, ext));
                    if ~isfile(dest_fqfn) && isfile(niis{i})
                        movefile(niis{i}, ipr.dest)
                    end
                end
            catch ME
                handwarning(ME)
            end
        end
        function j = nii2json(nii)
            assert(istext(nii))
            [pth,fp] = myfileparts(nii);
            j = fullfile(pth, strcat(fp, '.json'));
        end
        function g = original_filename_contains(str, opt)
            arguments
                str {mustBeText}
                opt.bids_folder {mustBeTextScalar} = 'rawdata'
                opt.modal_folder {mustBeTextScalar} = 'anat'
                opt.ext {mustBeTextScalar} = '.json'
                opt.ignorecase logical = false
            end

            p0 = pushd(fullfile(getenv('ADNI_HOME'), 'bids', opt.bids_folder));
            g = {};
            for g_ = glob(fullfile('sub-*', 'ses-*', opt.modal_folder, '*.json'))'
                j = fileread(g_{1});
                try
                    j = jsondecode(j);
                    if contains(j.ADNI_INFO.OriginalFilename, str, 'IgnoreCase', opt.ignorecase)
                        g__ = strrep(g_{1}, '.json', opt.ext);
                        g = [g g__]; %#ok<AGROW> 
                    end
                catch ME
                    handwarning(ME)
                    fprintf('filename:  %s\n', g__);
                end
            end
            popd(p0);
        end
        function sub = rid2subfold(r)
            persistent subfolds_            
            if isempty(subfolds_)
                tic
                subfolds_ = glob( ...
                    fullfile(getenv('ADNI_HOME'), 'bids', 'rawdata', 'sub-*'))';
                fprintf('rid2subfold:  '); toc
            end           
            assert(length(subfolds_) == 1660)
            sub = subfolds_(contains(subfolds_, r));
        end
        function tf = samedate(d1, d2)
            arguments
                d1 datetime
                d2 datetime
            end
            d1 = ensureTimeZone(d1);
            d2 = ensureTimeZone(d2);
            tf = duration(d1 - d2) < hours(24);
        end
        function rid = subfold2rid(s)
            re = regexp(s, 'sub-\d{3}S(?<rid>\d{4})', 'names');
            rid = str2double(re.rid);
        end             
        function m = viscode2months(carr)
            sarr = strrep(carr, "m", "");
            sarr = strrep(sarr, "bl", "0");
            for i = 1:length(sarr)
                m(i) = str2num(sarr(i)); %#ok<ST2NM,AGROW> 
            end
        end
    end

    %% PROTECTED

    properties (Access = protected)
        adni_demo_
        adni_merge_
    end

    methods (Access = protected)
        function this = AdniBids(varargin)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
