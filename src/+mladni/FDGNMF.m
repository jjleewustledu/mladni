classdef FDGNMF < handle & mladni.NMF
    %% 
    %  
    %  Created 23-Jun-2022 13:18:55 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.    
    
    methods
        function this = FDGNMF(varargin)
            this = this@mladni.NMF(varargin{:});            
        end
    end

    methods (Static)

        %% DEPRECATED
        %  Simple utilities for creating & adjustinng filenames, including BIDS filenames.

        function create_baseline_csv(fin, fout, varargin) 
            %% creates csvs of filenames for baseline* nmf datasets

            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'fin', @(x) isfile(x) && contains(x, '.csv'));
            addRequired(ip, 'fout', @(x) istext(x) && contains(x, '.csv'));
            addParameter(ip, 'modality', 'fdg', @istext);
            parse(ip, fin, fout, varargin{:});
            ipr = ip.Results;
            assert(~strcmp(ipr.fin, ipr.fout));
            ftmp = [tempname '.csv'];
            
            tbl_10 = readtable(ipr.fin, 'ReadVariableNames', false); % 10 vars split by "/"; var1 is nan
            tbl_10 = sortrows(tbl_10, 10);
            subs = {};
            fqfns = {};
            for row = 1:size(tbl_10,1)
                sub = tbl_10{row,7}{1};
                while ~contains(sub, subs)
                    subs = [subs; sub]; %#ok<*AGROW>
                    parts = tbl_10{row, 2:end}(:)';
                    fqfn = fullfile(filesep, parts{:});
                    fqfns = [fqfns; fqfn];
                end
            end            
            writetable(table(fqfns), ftmp, 'WriteVariableNames', false);
            
            % adjust csv for derivatives, proc, suffix
            mladni.FDGNMF.adjust_modality_csv(ftmp, ipr.fout, varargin{:});
            deleteExisting(ftmp)
        end
        function adjust_modality_csv(fin, fout, varargin)
            %% conveniently calls adjust_fdg_csv, adjust_t1w_csv.

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'modality', 'fdg', @istext);
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            switch lower(ipr.modality)
                case 'fdg'
                    mladni.FDGNMF.adjust_fdg_csv(fin, fout, varargin{:});
                case 't1w'
                    mladni.FDGNMF.adjust_t1w_csv(fin, fout, varargin{:});
                otherwise
                    error('mladni:ValueError', 'FDGNMF.adjust_modality_csv.ipr.modality: %s', ipr.modality);
            end
        end
        function adjust_fdg_csv(fin, fout, varargin)
            %% adjust BIDS names:  proc- and final suffix.

            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'fin', @(x) isfile(x) && contains(x, '.csv'));
            addRequired(ip, 'fout', @(x) istext(x) && contains(x, '.csv'))
            addParameter(ip, 'suffix', '_on_T1w_Warped_dlicv', @istext);
            addParameter(ip, 'proc', strcat(mladni.FDG.PROC, '-ponsvermis'), @istext);
            parse(ip, fin, fout, varargin{:});
            ipr = ip.Results;
            assert(~strcmp(ipr.fin, ipr.fout));
            
            tbl = readtable(ipr.fin, 'ReadVariableNames', false, 'Delimiter', ' ');
            c = cell(size(tbl));
            for row = 1:size(tbl,1)
                s = tbl{row,1};
                s = s{1};
                re = regexp(s, '\S+/sub-\d{3}S\d{4}_ses-\d{14}_\S*proc-(?<proc>[a-zA-Z\-]+)(?<tail>\S*).nii.gz', 'names');
                s = strrep(s, 'rawdata', 'derivatives');
                s = strrep(s, re.proc, ipr.proc);
                s = strrep(s, '.nii.gz', strcat(ipr.suffix, '.nii.gz'));
                c{row} = s;
            end
            writetable(table(c), ipr.fout, 'WriteVariableNames', false)
        end
        function adjust_t1w_csv(fin, fout, varargin)
            %% adjust BIDS names:  pve_idx

            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'fin', @(x) isfile(x) && contains(x, '.csv'));
            addRequired(ip, 'fout', @(x) istext(x) && contains(x, '.csv'));
            addParameter(ip, 'pve_idx', 1, @isscalar);
            parse(ip, fin, fout, varargin{:});
            ipr = ip.Results;
            assert(~strcmp(ipr.fin, ipr.fout));

            tree = readlines(fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'tree.log'), ...
                'EmptyLineRule', 'skip'); % string array
            tbl = readtable(ipr.fin, 'ReadVariableNames', false, 'Delimiter', ' ');
            c = {};
            for row = 1:size(tbl,1)
                s = tbl{row,1};
                s = s{1};
                pth = fileparts(s);
                pth = strrep(pth, 'rawdata', 'derivatives');
                re = regexp(s, '\S+/(?<sub>sub-\d{3}S\d{4})_(?<ses>ses-\d{14})_(?<tail>\S+).nii.gz', 'names');
                branch = tree(contains(tree, re.sub) & contains(tree, sprintf('_T1w_brain_pve_%i_Warped.nii.gz', ipr.pve_idx)));
                if isempty(branch)
                    continue
                end
                if length(branch) > 1
                    branch = branch(1);
                end
                fn = char(regexp(branch, 'sub-\d{3}S\d{4}\S+.nii.gz', 'match'));
                c = [c; fullfile(pth, fn)];
            end
            writetable(table(c), ipr.fout, 'WriteVariableNames', false);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
