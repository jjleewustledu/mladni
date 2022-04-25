classdef AdniBids
    %% provides supporting methods for class hierarchy
    %  
    %  Created 30-Dec-2021 16:40:32 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1809720 (R2021b) Update 1 for MACI64.  Copyright 2021 John J. Lee.
    
    methods (Static)
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
            %      bidsDir (folder):  default $SINGULARITY_HOME?ADNI/bids/rawdata

            import mladni.AdniBids.nii2json
            disp('Start mladni.AdniBids.batch_3dresample()')            
            t0 = tic;
            
            bidsDir = fullfile('/home/aris_data', 'ADNI_FDG', 'bids', 'rawdata', '');

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
            assert(contains(fn, '.nii'))
            if endsWith(fn, '.nii')
                fn_ = fn;
                fn = gzip(fn_);
                delete(fn_);
            end
            [pth,fp,e] = myfileparts(fn);
            re = regexp(fp, '(?<prefix>\S+)(?<suffix>(_T1\w*|_t1\w*|_pet))', 'names');
            fn1 = fullfile(pth, strcat(re.prefix, '_orient-rpi', re.suffix, e));
            cmd = sprintf('3dresample -debug 1 -orient rpi -prefix %s -input %s', fn1, fn);
            [~,r] = mlbash(cmd);            
            assert(isfile(fn1));
            
            % manage json
            j0 = fileread(strcat(myfileprefix(fn), '.json'));
            j1.afni_3dresample.cmd = cmd;
            j1.afni_3dresample.cmdout = r;
            jsonrecode(j0, j1, 'filenameNew', strcat(myfileprefix(fn1), '.json'));
            
            % clean previous
            deleteExisting(fn);
        end
        function dt = fqnifti2dt(fqniis)
            if ischar(fqniis) || isstring(fqniis)
                fqniis = {fqniis};
            end
            dt = NaT(size(fqniis));
            for fi = 1:length(fqniis)
                [~,s] = fileparts(fileparts(fileparts(fqniis{fi})));
                dt(fi) = mladni.AdniBids.ses2AcqDate(s); 
            end
        end
        function uid = fqnifti2uid(fqniis)
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

            import mladni.AdniBids.niis2datetimes
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
                dts = niis2datetimes(ipr.niis);
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
        function proc = nifti2proc(fn)
            %% char Description

            assert(istext(fn))
            re = regexp(fn, 'sub-\d{3}S\d{4}_ses-\d{14}_trc-\w+_proc-(?<proc>[CASUD]+)_pet', 'names');
            proc = re.proc;
        end
        function ses = nifti2ses(varargin)
            %% {8,14}-digit ses datetime

            ip = inputParser;
            addRequired(ip, 'fn', @istext)
            addOptional(ip, 'digits', 8, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;

            if 8 == ipr.digits
                re = regexp(ipr.fn, 'sub-\d{3}S\d{4}_(?<ses>ses-\d{8})\d{6}_\w+', 'names');
            else
                re = regexp(ipr.fn, 'sub-\d{3}S\d{4}_(?<ses>ses-\d{14})_\w+', 'names');
            end
            ses = re.ses;
        end
        function sub = nifti2sub(fn)
            assert(istext(fn))
            re = regexp(fn, '(?<sub>sub-\d{3}S\d{4})_\w+', 'names');
            sub = re.sub;
        end
        function dt = nii2datetime(nii)
            assert(istext(nii))
            re = regexp(nii, "sub-\d{3}S\d{4}_ses-(?<dt>\d{14})_\w+", "names");
            dt = datetime(re.dt, "InputFormat", "yyyyMMddHHmmss");
        end
        function dts = niis2datetimes(niis)
            assert(iscell(niis))
            for i = 1:length(niis)
                re = regexp(niis{i}, "sub-\d{3}S\d{4}_ses-(?<dt>\d{14})_\w+", "names");
                dts(i) = datetime(re.dt, "InputFormat", "yyyyMMddHHmmss"); %#ok<AGROW> 
            end
        end
        function j = nii2json(nii)
            assert(istext(nii))
            [pth,fp] = myfileparts(nii);
            j = fullfile(pth, strcat(fp, '.json'));
        end
        function d = proc2Description(p)
            assert(istext(p))
            switch upper(p)
                case 'CA'
                    d = 'Co-registered, Averaged';
                case 'CAS'
                    d = 'Coreg, Avg, Standardized Image and Voxel Size';
                case 'CASU'
                    d = 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution';
                case 'CD'
                    d = 'Co-registered Dynamic';
                otherwise
                    d = 'Unknown';
            end
        end
        function ad = ses2AcqDate(s)
            assert(istext(s))
            ad_ = regexp(s, 'ses-(?<y>\d{4})(?<m>\d{2})(?<d>\d{2})', 'names');
            ad = datetime(str2double(ad_.y), str2double(ad_.m), str2double(ad_.d));    
        end
        function ads = ses2AcqDatestr(s)
            assert(istext(s))
            ad_ = regexp(s, 'ses-(?<y>\d{4})(?<M>\d{2})(?<d>\d{2})(?<H>\d{2})(?<m>\d{2})(?<s>\d{2})', 'names');
            ads = sprintf('%s-%s-%s_%s_%s_%s.0', ad_.y, ad_.M, ad_.d, ad_.H, ad_.m, ad_.s);
        end
        function rid = sub2rid(s)
            re = regexp(s, 'sub-\d{3}S(?<rid>\d{4})', 'names');
            rid = str2double(re.rid);
        end
        function S = sub2Subject(s)
            assert(istext(s))
            s_ = regexp(s, 'sub-(?<n1>\d{3})S(?<n2>\d{4})', 'names');
            S = sprintf('%s_S_%s', s_.n1, s_.n2);
        end        
        function update_json_fdg(varargin)
            import mladni.AdniBids.sub2Subject
            import mladni.AdniBids.ses2AcqDate
            import mladni.AdniBids.proc2Description

            ip = inputParser;
            addRequired(ip, 'fn', @isfile)
            addRequired(ip, 't', @istable)
            addRequired(ip, 'sub', @istext)
            addRequired(ip, 'ses', @istext)
            addRequired(ip, 'proc', @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;

            s = sub2Subject(ipr.sub);
            ad = ses2AcqDate(ipr.ses);
            d = proc2Description(ipr.proc);
            t_ = ipr.t(contains(ipr.t.Subject, s) & ...
                       ipr.t.AcqDate == ad & ...
                       strcmp(ipr.t.Description, d), :);

            [pth,fp] = myfileparts(ipr.fn);
            fn = fullfile(pth, strcat(fp, '.json'));
            txt = jsonrecode(fn, table2struct(t_));
            try
                fid = fopen(fn, 'w');
                fprintf(fid, txt);
                fclose(fid);
            catch ME
                handwarning(ME)
                disp(fn)
            end
        end
        function t = update_table(~, varargin)
            import mladni.AdniBids.sub2Subject
            import mladni.AdniBids.ses2AcqDate
            import mladni.AdniBids.proc2Description
            import mladni.AdniBids.ses2AcqDatestr

            ip = inputParser;
            addRequired(ip, 'fn', @isfile)
            addRequired(ip, 't', @istable)
            addRequired(ip, 'sub', @istext)
            addRequired(ip, 'ses', @istext)
            addRequired(ip, 'proc', @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;            

            s = sub2Subject(ipr.sub);
            ad = ses2AcqDate(ipr.ses);
            d = proc2Description(ipr.proc);
            t = ipr.t;
            if ~contains(t.Properties.VariableNames, 'NativeFilename')                
                NativeFilename = cell(size(t.AcqDate));
                t = addvars(t, NativeFilename, 'Before', 'ImageDataID');
            end

            d1 = strrep(d, ' ', '_');
            ad1 = ses2AcqDatestr(ipr.ses);
            g = glob(fullfile(fileparts(ipr.fn), d1, ad1, 'I*', '*'));
            native_fn = fullfile(basename(pwd), g{1});

            t.NativeFilename(contains(t.Subject, s) & ...
                              t.AcqDate == ad & ...
                              strcmp(t.Description, d)) = {native_fn};
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

    methods (Access = protected)
        function this = AdniBids(varargin)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
