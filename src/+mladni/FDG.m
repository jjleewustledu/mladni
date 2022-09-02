classdef FDG < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% FDG uses object-oriented mlfsl.Flirt for registration and mlfourd.ImagingContext2 for actions on imaging data.  
    %  Most common use case:
    %  >> an_fdg = mlfourd.ImagingContext2(fq_filename)
    %  >> obj = mladni.FDG(an_fdg);
    %  >> obj.call()
    %  >> obj.view()
    %  
    %  Created 31-Dec-2021 00:02:37 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1809720 (R2021b) Update 1 for MACI64.  Copyright 2021 John J. Lee.

    methods (Static)
        function propcluster_tiny()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '6000'; % in MB; deepmrseg requires 10 GB; else 6 GB
            c.AdditionalProperties.Node = '';
            %c.AdditionalProperties.Partition = 'test';
            c.AdditionalProperties.WallTime = '4:00:00'; % 24 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function propcluster()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '6000'; % in MB; deepmrseg requires 10 GB; else 6 GB
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
        function [j,c] = parcluster_test()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.FDG.batch_test, 1, {31}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);            
        end
        function t = batch_test(len)
            t0 = tic;
            parfor idx = 1:len
                pause(10)
            end
            t = toc(t0);
        end
        
        function [j,c] = parcluster_tiny()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            %  --------------------------------------------------------------
            
            c = parcluster;
            disp(c.AdditionalProperties)

            globbing_file = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.mat');
            if ~isfile(globbing_file)
                j = c.batch(@mladni.FDG.batch_globbed, 1, {}, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);
                return
            end            
            ld = load(globbing_file);
            globbed = ld.globbed{1};
            globbed = globbed(1:3);
            j = c.batch(@mladni.FDG.batch_renorm_balanced, 1, {'globbed', globbed}, 'Pool', 3, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);      
        end
        function [j,c] = parcluster()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            %  --------------------------------------------------------------
            %  Use for call_resolved(), batch_revisit_pve1()
            
            c = parcluster;
            disp(c.AdditionalProperties)
            
            globbing_file = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.mat');
            if ~isfile(globbing_file)
                j = c.batch(@mladni.FDG.batch_globbed, 1, {}, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);
                return
            end            
            ld = load(globbing_file);
            for ji = 1:length(ld.globbed)
                j{ji} = c.batch(@mladni.FDG.batch_renorm_balanced, 1, {'globbed', ld.globbed{ji}}, 'Pool', 31, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);  %#ok<AGROW>
            end
        end
        function [j,c] = parcluster2()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            %  --------------------------------------------------------------
            %  Use for build_fdg_renormalized(), batch2()
            
            c = parcluster;
            disp(c.AdditionalProperties)
            
            globbing_file = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.mat');
            assert(isfile(globbing_file));
            ld = load(globbing_file);
            g = [ld.globbed{:}];
            g1 = cellfun(@(x) strrep(x, 'rawdata', 'derivatives'), ...
                g, 'UniformOutput', false);
            select = cellfun(@(x) isfile(x), g1);
            j = c.batch(@mladni.FDG.batch_renorm_balanced, 1, {'globbed', g(select)}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false); % {'len', []}, 'Pool', 31            
        end
        function parlocal2()
            %% renormalizes all FDG by {pons, cerebellar_vermis} per Susan Landau's methods

            globbing_file = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.mat');
            assert(isfile(globbing_file));
            ld = load(globbing_file);
            g = [ld.globbed{:}];
            g1 = cellfun(@(x) ...
                strrep(x, 'rawdata', 'derivatives'), ...
                g, 'UniformOutput', false);
            g2 = cellfun(@(x) ...
                strrep(x, '/home/aris_data/ADNI_FDG', '/mnt/CHPC_scratch/Singularity/ADNI'), ...
                g1, 'UniformOutput', false);
            select = cellfun(@(x) isfile(x), g2);
            globbed = g2(select);
            
            t0 = tic;
            parfor idx = 1:length(globbed)   
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});
                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.build_fdg_renormalized();
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            fprintf('mladni.FDG.parlocal2() completed in: ');
            toc(t0)
        end
        function parlocal3()
            %% binarize dlicv_detJ_Warped; use mask to extract brain from pet_on_T1w_Warped

            globbing_file = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.mat');
            assert(isfile(globbing_file));
            ld = load(globbing_file);
            g = [ld.globbed{:}];
            g1 = cellfun(@(x) ...
                strrep(x, 'rawdata', 'derivatives'), ...
                g, 'UniformOutput', false);
            g2 = cellfun(@(x) ...
                strrep(x, '/home/aris_data/ADNI_FDG', '/mnt/CHPC_scratch/Singularity/ADNI'), ...
                g1, 'UniformOutput', false);
            select = cellfun(@(x) isfile(x), g2);
            globbed = g2(select);
            
            t0 = tic;
            parfor idx = 1:length(globbed)   
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});
                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.build_pet_warped_brain();
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            fprintf('mladni.FDG.parlocal3() completed in: ');
            toc(t0)
        end
        function parlocal4()
            %% create images normalized by median values for cohort; enables fair use of FDG & PVE1 for NMF

            globbing_file = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.mat');
            assert(isfile(globbing_file));
            ld = load(globbing_file);
            g = [ld.globbed{:}];
            g1 = cellfun(@(x) ...
                strrep(x, 'rawdata', 'derivatives'), ...
                g, 'UniformOutput', false);
            g2 = cellfun(@(x) ...
                strrep(x, '/home/aris_data/ADNI_FDG', '/mnt/CHPC_scratch/Singularity/ADNI'), ...
                g1, 'UniformOutput', false);
            select = cellfun(@(x) isfile(x), g2);
            globbed = g2(select);
            
            workpth = '/mnt/CHPC_scratch/Singularity/ADNI/bids/derivatives';
            med = mlfourd.ImagingFormatContext2(fullfile(workpth, ...
                'all_trc-FDG_proc-CASU-ponsvermis-icv_orient-rpi_pet_on_T1w_Warped_dlicv_median.nii.gz'));
            med_fdg = median(med.img(med.img > eps));
            med = mlfourd.ImagingFormatContext2(fullfile(workpth, ...
                'all_acq-noaccel_proc-orig-n4_orient-rpi_T1w_brain_pve_1_detJ_Warped_median.nii.gz'));
            med_pve1 = median(med.img(med.img > eps));

            t0 = tic;
            parfor idx = 1:length(globbed)   
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});
                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.call_renorm_by_medians(med_fdg, med_pve1);
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            fprintf('mladni.FDG.parlocal4() completed in: ');
            toc(t0)
        end
        function t = batch_tiny(varargin)
            %% Requires completion of mladni.AdniBidsT1w.batch, mladniAdniBidsFdg.batch,
            %  and mladni.AdniBids.adni_3dresample to populate ADNI_FDG/bids/rawdata.
            %  Params:
            %      len (numeric):  # FDG sessions
            %      proc (text):  default 'CASU'
            
            ip = inputParser;
            addParameter(ip, 'len', 1:3, @isnumeric)
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            if isscalar(ipr.len)
                ipr.len = 1:ipr.len;
            end
            
            disp('Start mladni.FDG.batch_tiny()')        
            
            globbed = globT( ...
                fullfile( ...
                    '/home/aris_data/ADNI_FDG/bids/rawdata', ...
                    'sub-*', 'ses-*', 'pet', ...
                    sprintf('*-%s_*_pet.nii.gz', ipr.proc))); % *_proc-CASU_orient-rpi_pet
            globbed = globbed(ipr.len);
            %t = cell2table(globbed', 'VariableNames', {'rawdata_pet_filename'});
            %writetable(t, '/home/aris_data/ADNI_FDG/bids/derivatives/globbed_tiny.csv');
            fprintf('size(mladni.FDG.batch.globbed) -> %s\n', mat2str(size(globbed)));            
            
            t0 = tic;
            mladni.CHPC3.clean_tempdir();
            disp(ascol(globbed))
            parfor idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});
                    this = mladni.FDG(fdg_);
                    this.build_fdg_renormalized();
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end                     
            t = toc(t0);

            disp('mladni.FDG.batch_tiny() completed')
        end
        function globbed = batch_globbed(varargin)
            %% creates globbed.csv, globbed.mat
            
            ip = inputParser;
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            g = globT( ...
                fullfile( ...
                    '/home/aris_data/ADNI_FDG/bids/rawdata', ...
                    'sub-*', 'ses-*', 'pet', ...
                    sprintf('*-%s_*_pet.nii.gz', ipr.proc))); % *_proc-CASU_orient-rpi_pet
            len_part = ceil(length(g)/10);
            globbed = cell(1,10);
            for gi = 1:9
                globbed{gi} = g((gi-1)*len_part+1:gi*len_part);
            end
            globbed{10} = g(9*len_part+1:end);                
            
            t = cell2table(g', 'VariableNames', {'rawdata_pet_filename'});
            writetable( ...
                t, fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.csv'));
            save( ...
                fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.mat'), 'globbed');            
            fprintf('mladni.FDG.batch.globbed:\n')
            disp(globbed)
        end
        function t = batch(varargin)
            %% Requires completion of mladni.batch_globbed().
            %  Params:
            %      len (numeric):  # FDG sessions
            %      proc (text):  default 'CASU'
            
            ip = inputParser;
            addParameter(ip, 'globbed', {}, @iscell)
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mladni.FDG.batch()')        

            globbed = ipr.globbed;          
            %t = size(globbed); % DEBUG
            %disp(t) % DEBUG
            %return
            
            t0 = tic;
            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});
                    %derpth = strrep(fdg_.filepath, 'rawdata', 'derivatives');
                    %c{idx} = derpth; %% DEBUG
                    %continue

                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.call_resolve(); 
                        this.finalize(); 
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            t = toc(t0);

            disp('mladni.FDG.batch() completed')            
        end
        function t = batch2(varargin)
            %% Requires completion of mladni.batch_globbed().
            %  Params:
            %      len (numeric):  # FDG sessions
            %      proc (text):  default 'CASU'
            
            ip = inputParser;
            addParameter(ip, 'globbed', {}, @iscell)
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mladni.FDG.batch()')        

            globbed = ipr.globbed;          
            %t = size(globbed); % DEBUG
            %disp(t) % DEBUG
            %return
            
            t0 = tic;
            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});
                    %derpth = strrep(fdg_.filepath, 'rawdata', 'derivatives');
                    %c{idx} = derpth; %% DEBUG
                    %continue

                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.build_fdg_renormalized();
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            t = toc(t0);

            disp('mladni.FDG.batch() completed')            
        end
        function t = batch_renorm_balanced(varargin)
            %% Requires completion of mladni.batch_globbed().  Facilitates NMF for FDG & PVE1.
            %  Params:
            %      len (numeric):  # FDG sessions
            %      proc (text):  default 'CASU'
            
            ip = inputParser;
            addParameter(ip, 'globbed', {}, @iscell)
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mladni.FDG.batch_renorm_balanced()')        

            workpth = '/scratch/jjlee/Singularity/ADNI/bids/derivatives';
            med_fdg = mlfourd.ImagingContext2(fullfile(workpth, ...
                'all_trc-FDG_proc-CASU-ponsvermis_orient-rpi_pet_on_T1w_Warped_dlicv_median.nii.gz'));
            med_pve1 = mlfourd.ImagingContext2(fullfile(workpth, ...
                'all_acq-noaccel_proc-orig-n4_orient-rpi_T1w_brain_pve_1_Warped_median.nii.gz'));

            globbed = ipr.globbed;
            
            t0 = tic;
            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});

                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.call_renorm_balanced(med_fdg, med_pve1);
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            t = toc(t0);

            disp('mladni.FDG.batch_renorm_balanced() completed')  
        end
        function t = batch_renorm_by_medians(varargin)
            %% Requires completion of mladni.batch_globbed().
            %  Params:
            %      len (numeric):  # FDG sessions
            %      proc (text):  default 'CASU'
            
            ip = inputParser;
            addParameter(ip, 'globbed', {}, @iscell)
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mladni.FDG.batch_renorm_by_medians()')        

            workpth = '/scratch/jjlee/Singularity/ADNI/bids/derivatives';
            med = mlfourd.ImagingFormatContext2(fullfile(workpth, ...
                'all_trc-FDG_proc-CASU-ponsvermis-icv_orient-rpi_pet_on_T1w_Warped_dlicv_median.nii.gz'));
            med_fdg = median(med.img(med.img > eps));
            med = mlfourd.ImagingFormatContext2(fullfile(workpth, ...
                'all_acq-noaccel_proc-orig-n4_orient-rpi_T1w_brain_pve_1_detJ_Warped_median.nii.gz'));
            med_pve1 = median(med.img(med.img > eps));

            globbed = ipr.globbed;
            
            t0 = tic;
            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});

                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.call_renorm_by_medians(med_fdg, med_pve1);
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            t = toc(t0);

            disp('mladni.FDG.batch_renorm_by_medians() completed')  
        end
        function t = batch_revisit_pve1(varargin)
            %% Requires completion of mladni.batch_globbed().
            %  Params:
            %      len (numeric):  # FDG sessions
            %      proc (text):  default 'CASU'
            
            ip = inputParser;
            addParameter(ip, 'globbed', {}, @iscell)
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mladni.FDG.batch_revisit_pve1()')        

            globbed = ipr.globbed;
            
            t0 = tic;
            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});

                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.call_revisit_pve1();
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            t = toc(t0);

            disp('mladni.FDG.batch_revisit_pve1() completed')            
        end
        function t = batch_revisit_fdg(varargin)
            %% Requires completion of mladni.batch_globbed().
            %  Params:
            %      len (numeric):  # FDG sessions
            %      proc (text):  default 'CASU'
            
            ip = inputParser;
            addParameter(ip, 'globbed', {}, @iscell)
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mladni.FDG.batch_revisit_fdg()')        

            globbed = ipr.globbed;
            
            t0 = tic;
            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});

                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.call_revisit_fdg();
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            t = toc(t0);

            disp('mladni.FDG.batch_revisit_fdg() completed')    
        end

        function ic = add_proc(ic, proc0, proc1)
            assert(istext(proc0));
            assert(istext(proc1));
            if ~startsWith(proc1, '-')
                proc1 = strcat('-', proc1);
            end
            ic.fileprefix = strrep( ...
                ic.fileprefix, ...
                strcat('proc-', proc0), ...
                strcat('proc-', proc0, proc1));
        end
        function dt = AcqDate(fp)
            %% per mladni.AdniDemographics
            
            if isa(fp, 'mlio.IOInterface')
                fp = fp.fileprefix;
            end
            if contains(fp, filesep)
                [~,fp] = myfileparts(fp);
            end
            re = regexp(fp, 'sub-\d{3}S\d{4}_ses-(?<date>\d{8})', 'names');
            dt = datetime(re.date, 'InputFormat', 'yyyyMMdd');
        end        
        function ic = create_mask_for_nmf(varargin)

            ip = inputParser;
            addRequired(ip, 'fn_csv', @isfile);
            addParameter(ip, 'thresh', 0.5, @isscalar);
            parse(ip, varargin{:});
            ipr = ip.Results;

            filelist = readlines(ipr.fn_csv);

            ic = mlfourd.ImagingContext2(filelist{1});
            ic.selectNiftiTool();
            for fidx = 2:length(filelist)
                if isempty(filelist{fidx}); continue; end
                try
                    item = filelist{fidx};
                    if ~contains(item, 'orient-rpi_T1w_brain_pve_1_Warped.nii.gz')
                        pth = fileparts(item);
                        g = glob(fullfile(pth, '*orient-rpi_T1w_brain_pve_1_Warped.nii.gz'));
                        assert(~isempty(g))
                        item = g{1};
                    end
                    ic = ic + mlfourd.ImagingContext2(item);
                    ic.fileprefix = 'mladi_FDG_create_mask_for_nmf';
                catch ME
                    handwarning(ME);
                end
            end
            ic = ic/length(filelist);
            ic = ic.thresh(ipr.thresh);
            ic = ic.binarized();
            thr_ = strrep(num2str(ipr.thresh), '.', 'p');
            ic.fileprefix = sprintf('mladi_FDG_create_mask_for_nmf_thr%s', thr_);
            ic.save();
        end
        function m  = image_mass(ic)
            dV = voxelVolume(ic);
            ic1 = ic.thresh(0);
            sumDensities = dipsum(ic1);            
            m = sumDensities*dV;
        end
        function fn = json(obj)
            if isa(obj, 'mlfourd.ImagingContext2')
                fn = strcat(obj.fqfp, '.json');
                return
            end
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.json');
        end
        function jsonrecode(in, field, out)
            if istext(in)
                in = myfileprefix(in);
            end
            if isa(in, 'mlio.IOInterface')
                in = in.fqfp;
            end
            if istext(out)
                out = myfileprefix(out);
            end
            if isa(out, 'mlio.IOInterface')
                out = out.fqfp;
            end
            
            jsonrecode(strcat(in, '.json'), ...
                struct(clientname(true, 3), field), ...
                'filenameNew', strcat(out, '.json'));
        end
        function fn = mat(obj)
            if isa(obj, 'mlfourd.ImagingContext2')
                fn = strcat(obj.fqfp, '.mat');
                return
            end
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.mat');
        end
        function fn = niigz(obj)
            if isa(obj, 'mlfourd.ImagingContext2')
                fn = strcat(obj.fqfp, '.nii.gz');
                return
            end
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.nii.gz');
        end        
        function tbl = find_nans()
            %% find images with nans according to pons-vermis reference
            
            globbed = readlines(fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.csv'));
            ad = mladni.AdniDemographics();
            fdg1 = ad.table_fdg1();
            
            c = {};
            for gi = 1:length(globbed)
                re = regexp(globbed{gi}, ...
                    '/home/aris_data/ADNI_FDG/bids/rawdata/sub-(?<site>\d{3})S(?<rid>\d{4})/\S+', 'names');
                adniSub = sprintf('%s_S_%s', re.site, re.rid);
                select = contains(fdg1.Subject, adniSub);
                select1 = contains(fdg1.Description, 'Coreg, Avg, Standardized Image and Voxel Size');
                pv = fdg1.PonsVermis(select & select1);
                if isnan(pv)
                    stem = myfileprefix(globbed{gi});
                    stem = strrep(stem, 'rawdata', 'derivatives');
                    stem = strrep(stem, 'proc-CASU', 'proc-CASU-ponsvermis');
                    c = [c; glob([stem '*.nii.gz'])];   %#ok<AGROW>
                end
            end
            tbl = table(c);
        end
        function scrub_missing_from_csv(fn)
            assert(isfile(fn));
            assert(contains(hostname, 'cluster'));
            tbl = readtable(fn, 'ReadVariableNames', false, 'Delimiter', ' ');
            
            select = ~isfile(tbl.Var1);
            tbl(select,:) = [];
            if sum(select) > 0
                writetable(tbl, fn, 'WriteVariableNames', false);
            end
        end
        function scrub_nans_from_csv(fn)
            assert(isfile(fn));
            tbl = readtable(fn, 'ReadVariableNames', false, 'Delimiter', ' ');
            
            tbl_nans = readtable( ...
                fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed_nan.csv'), ...
                'ReadVariableNames', false, 'Delimiter', ' ');
            if contains(tbl.Var1{1}, 'derivatives')
                Var1 = cellfun(@(x) strrep(x, 'rawdata', 'derivatives'), tbl_nans.Var1, 'UniformOutput', false);
                tbl_nans = table(Var1);
            end
            if contains(tbl.Var1{1}, 'CASU-ponsvermis')
                Var1 = cellfun(@(x) strrep(x, 'CASU', 'CASU-ponsvermis'), tbl_nans.Var1, 'UniformOutput', false);
                tbl_nans = table(Var1);
            end
            if contains(tbl.Var1{1}, '_on_T1w_Warped.nii.gz')
                Var1 = cellfun(@(x) strrep(x, '.nii.gz', '_on_T1w_Warped.nii.gz'), tbl_nans.Var1, 'UniformOutput', false);
                tbl_nans = table(Var1);
            end
            
            select = contains(tbl.Var1, tbl_nans.Var1);
            tbl(select,:) = [];
            if sum(select) > 0
                writetable(tbl, fn, 'WriteVariableNames', false);
            end
        end
        function S  = Subject(fp)
            %% per mladni.AdniDemographics
            
            if isa(fp, 'mlio.IOInterface')
                fp = fp.fileprefix;
            end
            if contains(fp, filesep)
                [~,fp] = myfileparts(fp);
            end
            re = regexp(fp, 'sub-(?<pre>\d{3})S(?<rid>\d{4})_\S+', 'names');
            S = sprintf('%s_S_%s', re.pre, re.rid);
        end        
    end

    properties
        ants
        flirt_fdg
        flirt_t1w
    end

    properties (Dependent)
        blur % for co-registration of t1w with fdg
        debug % ~isempty(getenv('DEBUG'))
        Description

        atl
        atl_brain % skull-stripped brain only
        atl_mask 
        atl_mask1 
        
        fast_csf
        fast_csf_detJ_warped
        fast_csf_warped
        fast_gray
        fast_gray_detJ_warped
        fast_gray_warped
        fast_seg
        fast_white
        fast_white_detJ_warped        
        fast_white_warped
        
        fdg % mlfourd.ImagingContext2
        fdg_mask % mlfourd.ImagingContext2 
        fdg_mskt
        fdg_on_atl
        fdg_on_t1w % mlfourd.ImagingContext2         
        fdg_detJ_warped % warped to atl, weighted by det(J)
        fdg_detJ_warped_mask % generated from dlicv, binarized, 8mm blurring, thresh 0.1, binarized   
        fdg_proc % e.g., 'CASU' 
        fdg_reference % for renormalizing, e.g., ADNI, ponsvermis
        fdg_warped % warped to atl, not weighted by det(J)
        fdg_warped_dlicv % warped to atl, masked by dlicv, not weighted by det(J)
        
        t1w  
        t1w_blurred
        t1w_brain
        t1w_detJ
        t1w_dlicv
        t1w_dlicv_detJ_warped
        t1w_mskt
        t1w_n4
        t1w_on_atl
        t1w_warped
        t1w_0genericaffine
        t1w_1warp
        t1w_2warp % the composition of t1w_1warp \odot t1w_0genericaffine
    end

    methods

        %% GET, SET

        function g = get.blur(this)
            g = this.blur_;
        end
        function g = get.debug(~)
            g = ~isempty(getenv('DEBUG'));
        end
        function g = get.Description(this)
            switch (this.fdg_proc)
                case 'CA'
                    g = 'Co-registered, Averaged';
                case 'CD'
                    g = 'Co-registered Dynamic';
                case 'CAS'
                    g = 'Coreg, Avg, Standardized Image and Voxel Size';
                case 'CASU'
                    g = 'Coreg, Avg, Std Img and Vox Siz, Uniform Resolution';
                otherwise
                    error('mladni:ValueError', 'FDG.get.Description');
            end
        end
        
        function g = get.atl(this)
            g = this.atl_;
        end
        function g = get.atl_brain(this)
            fqfn = strcat(this.atl.fqfp, '_brain', '.nii.gz');
            g = mlfourd.ImagingContext2(fqfn);
        end
        function g = get.atl_mask(this)
            if ~isempty(this.atl_mask_)
                g = this.atl_mask_;
                return
            end
            this.atl_mask_ = mlfourd.ImagingContext2(strcat(this.atl.fqfp, '_brain_mask_dil', '.nii.gz'));
            g = this.atl_mask_;
        end
        function g = get.atl_mask1(this)
            if ~isempty(this.atl_mask1_)
                g = this.atl_mask1_;
                return
            end
            this.atl_mask1_ = mlfourd.ImagingContext2(strcat(this.atl.fqfp, '_brain_mask_dil1', '.nii.gz'));
            g = this.atl_mask1_;
        end
        
        function g = get.fast_csf(this)
            if ~isempty(this.fast_csf_)
                g = this.fast_csf_;
                return
            end
            this.fast_csf_ = mlfourd.ImagingContext2(strcat(this.t1w_brain.fqfp, '_pve_0', '.nii.gz'));
            g = this.fast_csf_;
        end
        function g = get.fast_csf_detJ_warped(this)
            if ~isempty(this.fast_csf_detJ_warped_)
                g = this.fast_csf_detJ_warped_;
                return
            end
            this.fast_csf_detJ_warped_ = mlfourd.ImagingContext2(strcat(this.fast_csf.fqfp, '_detJ_Warped.nii.gz'));
            g = this.fast_csf_detJ_warped_;
        end
        function g = get.fast_csf_warped(this)
            if ~isempty(this.fast_csf_warped_)
                g = this.fast_csf_warped_;
                return
            end
            this.fast_csf_warped_ = mlfourd.ImagingContext2(strcat(this.fast_csf.fqfp, '_Warped.nii.gz'));
            g = this.fast_csf_warped_;
        end
        function g = get.fast_gray(this)
            if ~isempty(this.fast_gray_)
                g = this.fast_gray_;
                return
            end
            this.fast_gray_ = mlfourd.ImagingContext2(strcat(this.t1w_brain.fqfp, '_pve_1', '.nii.gz'));
            g = this.fast_gray_;
        end
        function g = get.fast_gray_detJ_warped(this)
            if ~isempty(this.fast_gray_detJ_warped_)
                g = this.fast_gray_detJ_warped_;
                return
            end
            this.fast_gray_detJ_warped_ = mlfourd.ImagingContext2(strcat(this.fast_gray.fqfp, '_detJ_Warped.nii.gz'));
            g = this.fast_gray_detJ_warped_;
        end
        function g = get.fast_gray_warped(this)
            if ~isempty(this.fast_gray_warped_)
                g = this.fast_gray_warped_;
                return
            end
            this.fast_gray_warped_ = mlfourd.ImagingContext2(strcat(this.fast_gray.fqfp, '_Warped.nii.gz'));
            g = this.fast_gray_warped_;
        end
        function g = get.fast_seg(this)
            if ~isempty(this.fast_seg_)
                g = this.fast_seg_;
                return
            end
            this.fast_seg_ = mlfourd.ImagingContext2(strcat(this.t1w_brain.fqfp, '_seg', '.nii.gz'));
            g = this.fast_seg_;
        end
        function g = get.fast_white(this)
            if ~isempty(this.fast_white_)
                g = this.fast_white_;
                return
            end
            this.fast_white_ = mlfourd.ImagingContext2(strcat(this.t1w_brain.fqfp, '_pve_2', '.nii.gz'));
            g = this.fast_white_;
        end
        function g = get.fast_white_detJ_warped(this)
            if ~isempty(this.fast_white_detJ_warped_)
                g = this.fast_white_detJ_warped_;
                return
            end
            this.fast_white_detJ_warped_ = mlfourd.ImagingContext2(strcat(this.fast_white.fqfp, '_detJ_Warped.nii.gz'));
            g = this.fast_white_detJ_warped_;
        end
        function g = get.fast_white_warped(this)
            if ~isempty(this.fast_white_warped_)
                g = this.fast_white_warped_;
                return
            end
            this.fast_white_warped_ = mlfourd.ImagingContext2(strcat(this.fast_white.fqfp, '_Warped.nii.gz'));
            g = this.fast_white_warped_;
        end
  
        function g = get.fdg(this)
            g = this.fdg_;
        end
        function g = get.fdg_mask(this)
            if ~isempty(this.fdg_mskt_)
                g = this.fdg_mskt_;
                return
            end
            pwd0 = pushd(this.fdg.filepath);
            this.fdg_mskt_ = this.fdg.thresh(0.0001*dipmax(this.fdg));
            this.fdg_mskt_ = this.fdg_mskt_.binarized();
            this.fdg_mskt_.fileprefix = strcat(this.fdg.fileprefix, '_msk');
            this.fdg_mskt_.save();
            g = this.fdg_mskt_;
            popd(pwd0);
        end
        function g = get.fdg_mskt(this)
            if ~isempty(this.fdg_mskt_)
                g = this.fdg_mskt_;
                return
            end
            pwd0 = pushd(this.fdg.filepath);
            this.fdg_mskt_ = mlfourd.ImagingContext2( ...
                mlfsl.Flirt.msktgen(this.fdg.fqfn, 'cost', 'mutualinfo', 'dof', 6));
            g = this.fdg_mskt_;
            popd(pwd0);
        end
        function g = get.fdg_on_atl(this)
            if ~isempty(this.fdg_on_atl_)
                g = this.fdg_on_atl_;
                return
            end
            fqfp = this.in_derivatives_path(this.fdg.fqfp);
            this.fdg_on_atl_ = mlfourd.ImagingContext2(strcat(fqfp, '_on_atl.nii.gz'));
            g = this.fdg_on_atl_;
        end
        function g = get.fdg_on_t1w(this)
            if ~isempty(this.fdg_on_t1w_)
                g = this.fdg_on_t1w_;
                return
            end
            fqfp = this.in_derivatives_path(this.fdg.fqfp);
            this.fdg_on_t1w_ = mlfourd.ImagingContext2(strcat(fqfp, '_on_T1w.nii.gz'));
            if strcmp(this.fdg_reference, 'ponsvermis')
                this.fdg_on_t1w_ = this.add_proc(this.fdg_on_t1w_, this.fdg_proc, 'ponsvermis');
            end
            if strcmp(this.fdg_reference, 'ponsvermis-icv')
                this.fdg_on_t1w_ = this.add_proc(this.fdg_on_t1w_, this.fdg_proc, 'ponsvermis-icv');
            end
            g = this.fdg_on_t1w_;
        end
        function g = get.fdg_detJ_warped(this)
            if ~isempty(this.fdg_detJ_warped_)
                g = this.fdg_detJ_warped_;
                return
            end
            fqfp = this.in_derivatives_path(this.fdg_on_t1w.fqfp);
            this.fdg_detJ_warped_ = mlfourd.ImagingContext2(strcat(fqfp, '_detJ_Warped.nii.gz'));
            g = this.fdg_detJ_warped_;
        end        
        function g = get.fdg_detJ_warped_mask(this)
            if ~isempty(this.fdg_detJ_warped_mask_)
                g = this.fdg_detJ_warped_mask_;
                return
            end
            fqfp = this.fdg_detJ_warped.fqfp;
            this.fdg_detJ_warped_mask_ = mlfourd.ImagingContext2(strcat(fqfp, '_mask.nii.gz'));
            g = this.fdg_detJ_warped_mask_;
        end  
        function g = get.fdg_proc(this)
            g = this.fdg_proc_;
        end
        function g = get.fdg_reference(this)
            g = this.fdg_reference_;
        end
        function     set.fdg_reference(this, s)
            assert(contains(s, {'ADNI', 'ponsvermis', 'ponsvermis-icv'}, 'IgnoreCase', false));
            this.fdg_reference_ = s;

            % reset caches
            this.fdg_on_t1w_ = [];
            this.fdg_detJ_warped_ = [];
            this.fdg_warped_ = [];
        end
        function g = get.fdg_warped(this)
            if ~isempty(this.fdg_warped_)
                g = this.fdg_warped_;
                return
            end
            fqfp = this.in_derivatives_path(this.fdg_on_t1w.fqfp);
            this.fdg_warped_ = mlfourd.ImagingContext2(strcat(fqfp, '_Warped.nii.gz'));
            g = this.fdg_warped_;
        end 
        function g = get.fdg_warped_dlicv(this)
            if ~isempty(this.fdg_warped_dlicv_)
                g = this.fdg_warped_dlicv_;
                return
            end
            fqfp = this.in_derivatives_path(this.fdg_on_t1w.fqfp);
            this.fdg_warped_dlicv_ = mlfourd.ImagingContext2(strcat(fqfp, '_Warped_dlicv.nii.gz'));
            g = this.fdg_warped_dlicv_;
        end
        
        function g = get.t1w(this)
            g = this.t1w_;
        end
        function g = get.t1w_blurred(this)
            if ~isempty(this.t1w_blurred_)
                g = this.t1w_blurred_;
                return
            end
            if ~isempty(this.t1w_n4_)
                this.t1w_blurred_ = this.t1w_n4.blurred(this.blur);
            else
                this.t1w_blurred_ = this.t1w.blurred(this.blur);
            end
            g = this.t1w_blurred_;
            if ~isfile(g.fqfn)
                g.save();
            end
        end
        function g = get.t1w_brain(this)
            if ~isempty(this.t1w_brain_)
                g = this.t1w_brain_;
                return
            end
            if isfile(this.t1w_n4.fqfn) && isfile(this.t1w_dlicv.fqfn)
                fqfp = this.t1w_n4.fqfp;
                this.t1w_brain_ = this.t1w_n4 .* this.t1w_dlicv;
            else                
                fqfp = this.in_derivatives_path(this.t1w.fqfp);
                this.t1w_brain_ = this.t1w .* this.t1w_mskt;
            end
            this.t1w_brain_.fqfp = strcat(fqfp, '_brain.nii.gz');
            this.t1w_brain_.save();
            g = this.t1w_brain_;
        end
        function g = get.t1w_detJ(this)
            if ~isempty(this.t1w_detJ_)
                g = this.t1w_detJ_;
                return
            end
            fqfp = this.t1w_brain.fqfp;
            this.t1w_detJ_ = mlfourd.ImagingContext2(strcat(fqfp, '_detJ.nii.gz'));
            g = this.t1w_detJ_;
        end
        function g = get.t1w_dlicv(this)
            if ~isempty(this.t1w_dlicv_)
                g = this.t1w_dlicv_;
                return
            end
            fqfp = this.t1w_n4.fqfp;
            this.t1w_dlicv_ = mlfourd.ImagingContext2(strcat(fqfp, '_dlicv', '.nii.gz'));
            g = this.t1w_dlicv_;
        end
        function g = get.t1w_dlicv_detJ_warped(this)
            if ~isempty(this.t1w_dlicv_detJ_warped_)
                g = this.t1w_dlicv_detJ_warped_;
                return
            end
            fqfp = this.t1w_dlicv.fqfp;
            this.t1w_dlicv_detJ_warped_ = mlfourd.ImagingContext2(strcat(fqfp, '_detJ_Warped.nii.gz'));
            g = this.t1w_dlicv_detJ_warped_;
        end
        function g = get.t1w_mskt(this)
            if ~isempty(this.t1w_mask_)
                g = this.t1w_mask_;
                return
            end   
            pwd0 = pushd(this.t1w.filepath);         
            this.t1w_mask_ = mlfourd.ImagingContext2(mlfsl.Flirt.msktgen(this.t1w.fqfn));
            g = this.t1w_mask_;
            popd(pwd0);
        end
        function g = get.t1w_n4(this)
            if ~isempty(this.t1w_n4_)
                g = this.t1w_n4_;
                return
            end
            fqfp = this.t1w.fqfp;
            %if ~contains(fqfp, 'n3') && ~contains(fqfp, 'n4')
            this.t1w_n4_ = mlfourd.ImagingContext2( ...
                strcat(strrep(fqfp, '_orient', '-n4_orient'), '.nii.gz'));
            g = this.t1w_n4_;
            %else
            %    g = this.t1w;
            %end
        end        
        function g = get.t1w_on_atl(this)
            if ~isempty(this.t1w_on_atl_)
                g = this.t1w_on_atl_;
                return
            end
            fqfp = this.in_derivatives_path(this.t1w.fqfp);
            this.t1w_on_atl_ = mlfourd.ImagingContext2(strcat(fqfp, '_on_atl.nii.gz'));
            g = this.t1w_on_atl_;
        end
        function g = get.t1w_warped(this)
            if ~isempty(this.t1w_warped_)
                g = this.t1w_warped_;
                return
            end
            fqfp = this.t1w_brain.fqfp;
            this.t1w_warped_ = mlfourd.ImagingContext2(strcat(fqfp, '_Warped.nii.gz'));
            g = this.t1w_warped_;
        end
        function g = get.t1w_0genericaffine(this)
            if ~isempty(this.t1w_0genericaffine_)
                g = this.t1w_0genericaffine_;
                return
            end
            fqfp = this.t1w_brain.fqfp;
            this.t1w_0genericaffine_ = mlfourd.ImagingContext2(strcat(fqfp, '_0GenericAffine.nii.gz'));
            g = this.t1w_0genericaffine_;
        end
        function g = get.t1w_1warp(this)
            if ~isempty(this.t1w_1warp_)
                g = this.t1w_1warp_;
                return
            end
            fqfp = this.t1w_brain.fqfp;
            this.t1w_1warp_ = mlfourd.ImagingContext2(strcat(fqfp, '_1Warp.nii.gz'));
            g = this.t1w_1warp_;
        end
        function g = get.t1w_2warp(this)
            %% the composition of t1w_1warp \odot t1w_0genericaffine
            
            if ~isempty(this.t1w_2warp_)
                g = this.t1w_2warp_;
                return
            end
            fqfp = this.t1w_brain.fqfp;
            this.t1w_2warp_ = mlfourd.ImagingContext2(strcat(fqfp, '_2Warp.nii.gz'));
            g = this.t1w_2warp_;
        end

        %%

        function out  = antsApplyWarpToAtl(this, in)
            if ~isfile(in.fqfn)
                in.save();
            end
            out = mlfourd.ImagingContext2( ...
                this.ants.antsApplyTransforms(in, this.atl_brain, this.t1w_brain));
            out.selectNiftiTool();
        end
        function out  = antsApplyWarpToAtlNearNeigh(this, in)
            if ~isfile(in.fqfn)
                in.save();
            end
            out = mlfourd.ImagingContext2( ...
                this.ants.antsApplyTransforms(in, this.atl_brain, this.t1w_brain, ...
                '-n NearestNeighbor'));
            out.selectNiftiTool();
        end
        function this = applyXfm_fdg2atl(this)
            %% #APPLYXFM_FDG2ATL

            assert(~isempty(this.flirt_t1w));
            assert(isfile(this.niigz(this.fdg_on_t1w)));
            
            flirt_on_atl_ = copy(this.flirt_t1w);
            flirt_on_atl_.in = this.fdg_on_t1w;
            flirt_on_atl_.ref = this.atl;
            flirt_on_atl_.out = this.niigz(this.fdg_on_atl);
            flirt_on_atl_.applyXfm();

            if isfile(strcat(this.fdg_on_t1w.fqfp, ".json"))
                j0 = fileread(strcat(this.fdg_on_t1w.fqfp, ".json"));
                [~,j1] = this.flirt_t1w.cost_final();
                jsonrecode(j0, j1, 'filenameNew', this.json(flirt_on_atl_.out));
            end
        end
        function out  = applyXfm_fast2atl(this, in)
            %% #APPLYXFM_FAST2ATL
            %  Args:
            %      in (any): any imaging object in the space of this.t1w.

            assert(isfile(this.niigz(in)));            
            
            in = mlfourd.ImagingContext2(in);
            out = mlfourd.ImagingContext2(strcat(in.fqfp, '_on_atl.nii.gz'));
            if ~isempty(this.flirt_t1w)
                flirt_on_atl_ = copy(this.flirt_t1w);
                flirt_on_atl_.in = in;
                flirt_on_atl_.ref = this.atl;
                flirt_on_atl_.out = this.fast_gray;
            else
                flirt_on_atl_ = mlfsl.Flirt( ...
                    'in', in, ...
                    'ref', this.atl, ...
                    'out', this.fast_gray, ...
                    'omat', this.mat(this.t1w_on_atl), ...
                    'interp', 'trilinear');
            end
            flirt_on_atl_.applyXfm();
        end
        function out  = applyXfm_t1wlike2atl(this, in)
            %% #APPLYXFM_T1WLIKE2ATL
            %  Args:
            %      in (any): any imaging object in the space of this.t1w.

            assert(~isempty(this.flirt_t1w));
            assert(isfile(this.niigz(in)));            
            
            in = mlfourd.ImagingContext2(in);
            out = mlfourd.ImagingContext2(strcat(in.fqfp, '_on_atl.nii.gz'));
            flirt_on_atl_ = copy(this.flirt_t1w);
            flirt_on_atl_.in = in;
            flirt_on_atl_.ref = this.atl;
            flirt_on_atl_.out = this.niigz(out);
            flirt_on_atl_.applyXfm();

            if isfile(strcat(in.fqfp, ".json"))
                j0 = fileread(strcat(in.fqfp, ".json"));
                [~,j1] = this.flirt_t1w.cost_final();
                jsonrecode(j0, j1, 'filenameNew', this.json(flirt_on_atl_.out));
            end
        end        
        function this = build_deepmrseg_warped(this) 
            % warp to atlas    

            t1w_dlicv_w = mlfourd.ImagingContext2( ...
                this.ants.antsApplyTransforms(this.t1w_dlicv, this.atl_brain, this.t1w_brain));            
            this.t1w_dlicv_detJ_warped_ = this.t1w_detJ .* t1w_dlicv_w;
            this.t1w_dlicv_detJ_warped_.fqfp = strcat(this.t1w_dlicv.fqfp, '_detJ_Warped');
            this.t1w_dlicv_detJ_warped_.save();
            
            mladni.FDG.jsonrecode( ...
                t1w_dlicv_w, ...
                struct('state_changes', 'this.t1w_dlicv_detJ_warped_ = this.t1w_detJ .* t1w_dlicv_w', ...
                       'image_mass', this.image_mass(this.t1w_dlicv_detJ_warped_)), ...
                this.t1w_dlicv_detJ_warped_);

            deleteExisting(t1w_dlicv_w);
        end
        function this = build_fast_warped2(this)
            % warp all segs to atlas    
            
            for pve = {'fast_csf', 'fast_gray', 'fast_white'}                
                this.([pve{1} '_warped_']) = this.build_fast_warped_seg2(this.(pve{1}));
            end
        end
        function inw  = build_fast_warped_seg2(this, in)
            if ~isfile(in.fqfn)
                in.save();
            end
            inw = mlfourd.ImagingContext2( ...
                this.ants.antsApplyTransforms(in, this.atl_brain, this.t1w_brain));

            mladni.FDG.jsonrecode( ...
                inw, ...
                struct('state_changes', 'no use of this.t1w_detJ', ...
                       'image_mass', this.image_mass(inw)), ...
                inw);
        end
        function this = build_fdg_renormalized(this)
            %% renormalizes all FDG by {pons, cerebellar_vermis} per Susan Landau's methods
            
            try
                for ic = {this.fdg_warped} % this.fdg_on_t1w, this.fdg_detJ_warped
                    suvr = this.suvr_pons_vermis();
                    assert(suvr > 0.5 && suvr < 2);
                    icv = this.icv_from_imaging(this.t1w_dlicv);
                    
                    ic_ = ic{1};
                    ic_.selectNiftiTool();
                    ic_ = ic_./suvr;
                    ic_ = ic_./icv;
                    ic_ = ic_ .* binarized(this.t1w_dlicv_detJ_warped);
                    ic_.fileprefix = strrep(ic{1}.fileprefix, ...
                        strcat('proc-', this.fdg_proc), ...
                        strcat('proc-', this.fdg_proc, '-ponsvermis-icv'));
                    ic_.fileprefix = strcat(ic_.fileprefix, '_dlicv');
                    ic_.save();
                    
                    mladni.FDG.jsonrecode( ...
                        ic{1}, ...
                        struct('suvr_pons_vermis', suvr, 'icv', icv), ...
                        ic_);
                end
            catch ME % suvr not always available
                handwarning(ME)
            end
        end
        function this = build_fdg_warped(this)
            fdgw = mlfourd.ImagingContext2( ...
                this.ants.antsApplyTransforms(this.fdg_on_t1w, this.atl_brain, this.t1w_brain));
            fdgw = fdgw.thresh(0);
            this.fdg_detJ_warped_ = this.t1w_detJ .* fdgw;
            this.fdg_detJ_warped_.fqfp = strcat(this.fdg_on_t1w.fqfp, '_detJ_Warped');
            this.fdg_detJ_warped_.save();
            
            mladni.FDG.jsonrecode(...
                fdgw, ...
                struct('state_changes', 'fdgw = fdgw.thresh(0); this.fdg_detJ_warped_ = this.t1w_detJ .* fdgw', ...
                       'image_activity', this.image_mass(this.fdg_detJ_warped_)), ...
                this.fdg_detJ_warped_);

            deleteExisting(fdgw);
        end
        function this = build_nmf_inputs(this)
            %% builds fdg_detJ_warped_mask
            
            msk = this.t1w_dlicv_detJ_warped;
            msk = msk.binarized();
            msk.fqfp = strcat(this.fdg_detJ_warped.fqfp, '_mask'); 
            msk.save();
            msk = msk.blurred(8);
            
            for thr = 0.5:-0.1:0.1
                msk_ = msk.thresh(thr);
                msk_ = msk_.binarized();
                msk_.fqfp = strcat(this.fdg_detJ_warped.fqfp, sprintf('_mask_%g', thr));
                msk_.save();
            end  
            this.fdg_detJ_warped_mask_ = msk_;
        end
        function this = build_pet_warped_brain(this)
            %% binarize dlicv_detJ_Warped; use mask to extract brain from pet_on_T1w_Warped

            try
                bin = this.t1w_dlicv_detJ_warped.binarized();
                fdg_fqfn0 = this.fdg_warped.fqfn;
                this.fdg_reference = 'ponsvermis';
                if ~isfile(this.fdg_warped)
                    return
                end
                brain = this.fdg_warped .* bin;
                brain.fqfp = strcat(this.fdg_warped.fqfp, '_dlicv');
                brain.save();

                mladni.FDG.jsonrecode( ...
                    fdg_fqfn0, ...
                    struct('build_pet_warped_brain', 'bin = this.t1w_dlicv_detJ_warped.binarized()'), ...
                    brain);
            catch ME
                handwarning(ME)
            end

        end
        function this = call_renorm_balanced(this, med_fdg, med_pve1)
            %% extends call_renorm_by_medians by also weighting fdg intensities by # fdg voxels

            assert(isa(med_fdg, 'mlfourd.ImagingContext2'));
            assert(isa(med_pve1, 'mlfourd.ImagingContext2'));
            med_fdg_scalar = double(med_fdg.maskedMaths(med_fdg.binarized(), @median)); % median of img by binary mask
            med_pve1_scalar = double(med_pve1.maskedMaths(med_pve1.binarized(), @median));
            this.fdg_reference = 'ponsvermis';

            % FDG
            med_fdg_bin = med_fdg.binarized(); % # voxels
            med_pve1_bin = med_pve1.binarized(); % # voxels
            norm_fdg = dipsum(med_fdg_bin)/dipsum(med_pve1_bin)*med_fdg_scalar;
            ic_ = this.fdg_warped_dlicv;
            ic = ic_./norm_fdg;
            ic.fileprefix = strcat(ic_.fileprefix, strrep(sprintf('_renorm%g', norm_fdg), '.', 'p'));
            ic.save();
            mladni.FDG.jsonrecode( ...
                this.fdg_warped_dlicv, ...
                struct('state_changes', sprintf('renorm by dipsum(med_fdg_bin)/dipsum(med_pve1_bin)*med_fdg_scalar->%g', norm_fdg)), ...
                ic);

            % PVE1
            ic_ = this.fast_gray_detJ_warped;
            ic = ic_./med_pve1_scalar;
            ic.fileprefix = strcat(ic_.fileprefix, strrep(sprintf('_renorm%g', med_pve1_scalar), '.', 'p'));
            ic.save();
            mladni.FDG.jsonrecode( ...
                this.fast_gray_detJ_warped, ...
                struct('state_changes', sprintf('renorm by pve1 median->%g', med_pve1_scalar)), ...
                ic);
        end
        function this = call_renorm_by_medians(this, med_fdg, med_pve1)
            assert(isscalar(med_fdg));
            assert(isscalar(med_pve1));
            this.fdg_reference = 'ponsvermis';

            % FDG
            ic_ = this.fdg_warped_dlicv;
            ic = ic_./med_fdg;
            ic.fileprefix = strcat(ic_.fileprefix, strrep(sprintf('_renorm%g', med_fdg), '.', 'p'));
            ic.save();
            mladni.FDG.jsonrecode( ...
                this.fdg_warped_dlicv, ...
                struct('process', sprintf('renorm by median->%g', med_fdg)), ...
                ic);

            % PVE1
            ic_ = this.fast_gray_detJ_warped;
            ic = ic_./med_pve1;
            ic.fileprefix = strcat(ic_.fileprefix, strrep(sprintf('_renorm%g', med_pve1), '.', 'p'));
            ic.save();
            mladni.FDG.jsonrecode( ...
                this.fast_gray_detJ_warped, ...
                struct('process', sprintf('renorm by median->%g', med_pve1)), ...
                ic);
        end
        function this = call_resolve(this)
            %% #CALL_RESOLVE
            
            this = this.N4BiasFieldCorrection();
            this = this.deepmrseg_apply();
            this = this.resolve_fdg2t1w();
            this = this.fast();
            this = this.CreateJacobianDeterminantImages();
            this = this.build_deepmrseg_warped();
            %this = this.build_fast_warped();
            this = this.build_fdg_warped();
            %this = this.build_nmf_inputs();
            this = this.build_fdg_renormalized();
            %this = this.save_qc();
        end
        function this = call_revisit_fdg(this)
            %  - renormalize fdg by pons-vermis & intracranial volume

            this = this.build_fdg_renormalized();
        end
        function this = call_revisit_pve1(this)
            %% #CALL_REVISIT_PVE1
            %
            %  Aris to me, July 28th.
            %  OK. some additional thoughts from a discussion with Janine and Petra regarding smoothing and thresholding 
            %  the pve estimates. What Janine mentioned, which is reasonable, is that smoothing can create two problems 
            %      1) spill out signal as we see happening in the wm; but also 
            %      2) "corrupting" gm signal by "smoothing in" signal from background (i.e., zero values). 
            %  Thus, it might be better to smooth modulated pve values (in which case 0 background values will have 
            %  smaller influence) and used thresholded (non-modulated) pve values to create the mask needed as input for
            %  the nmf. The latter will be more strict and will allow us to focus on the gm ribbon avoiding some of the
            %  wm regions where we see loadings (e.g., component #7).            

            for t = {'gray'} % 'csf', 'white'

                %% make mask for command-line NMF, generating [this.fast_gray_Warped.nii.gz']

                % warp pve for externally generating group mask, preferably with subsequent thresh, binarization
                tiss__ = this.(strcat('fast_', t{1}));
                this.antsApplyWarpToAtl(tiss__);

                %% GM input to NMF is [this.fast_gray_detJ_Warped.nii.gz'], renormed by ICV in L

                % renorm pve by (intracranial volume from fast seg)/1e6
                icv_ = this.fast_seg_binary_for('icv');
                norm_ = dipsum(icv_)*prod(icv_.imagingFormat.mmppix)/1e6; % liters
                tiss_ = this.(strcat('fast_', t{1}, '_detJ_warped'));
                tiss_ = tiss_ ./ norm_;
                
                % smooth
                tiss_ = tiss_.blurred(7.9373); % sigma(8)^2 - sigma(1)^2 = sigma(7.9373)^2
                
                % save
                lbl = sprintf('fast_%s_detJ_warped', t{1});
                lbl_ = strcat(lbl, '_');
                tiss_.fileprefix = this.(lbl).fileprefix;
                tiss_.save();
                mladni.FDG.jsonrecode( ...
                    this.t1w, ...
                    struct('state_changes', sprintf( ...
                           'norm_ ~ icv_; blurred(1/norm_ .* antsApplyWarpToAtl(this, this.fast_%s_detJ_warped), 7.9373)', t{1}), ...
                           'image_mass', this.image_mass(tiss_)), ...
                    tiss_);
                this.(lbl_) = tiss_;
            end
        end
        function this = CreateJacobianDeterminantImages(this)
            %% #CREATEJACOBIANDETERMINANTIMAGE creates warping files, t1w_detJ.
            
            this.t1w_warped_ = mlfourd.ImagingContext2( ...
                this.ants.antsRegistrationSyNQuick(this.atl_brain, this.t1w_brain));
            this.ants.antsApplyTransforms2(this.t1w_brain, this.atl_brain, this.t1w_brain);
            this.ants.CreateJacobianDeterminantImage(this.t1w_2warp, this.t1w_detJ);
        end
        function this = deepmrseg_apply(this)
            if isfile(this.t1w_dlicv.fqfn)
                return
            end
            
            sif = fullfile(getenv('SINGULARITY_HOME'), 'deepmrseg_image_20220515.sif');
            cmd = sprintf('singularity exec --bind %s:/data %s "deepmrseg_apply" "--task" "dlicv" "--inImg" "/data/%s" "--outImg" "/data/%s"', ...
                this.t1w_n4.filepath, sif, this.t1w_n4.filename, this.t1w_dlicv.filename);
            mlbash(cmd);

            mladni.FDG.jsonrecode( ...
                this.t1w_n4, ...
                struct('bash', cmd, ...
                       'image_mass', this.image_mass(this.t1w_dlicv)), ...
                this.t1w_dlicv);
        end
        function fn   = duplicate_orig_file(~, obj)
            if ~istext(obj)
                obj = obj.fqfilename;
            end
            fn = strrep(obj, 'rawdata', 'derivatives');
            fn = strrep(fn, '_orient-rpi-std', '_orient-rpi');
        end
        function this = fast(this)
            pwd0 = pushd(this.t1w_brain.filepath);
            
            cmd = sprintf('fast -o %s %s', this.t1w_brain.fqfp, this.t1w_brain.fqfn);
            mlbash(cmd);
            
            for pve = {'fast_csf', 'fast_gray', 'fast_white'}
                mladni.FDG.jsonrecode( ...
                    this.t1w_brain, ...
                    struct('bash', cmd, ...
                           'image_mass', this.image_mass(this.(pve{1}))), ...
                    this.(pve{1}));
            end
            
            popd(pwd0);
        end
        function ic   = fast_pve_for(this, tissue)
            switch tissue
                case {'csf', 0}
                    ic = copy(this.fast_csf);
                case {'gray', 1}
                    ic = copy(this.fast_gray);
                case {'white', 2}
                    ic = copy(this.fast_white);
                otherwise
                    error('mladni:ValueError', 'FDG.fast_pve_for()')
            end
        end
        function ic   = fast_seg_binary_for(this, tissue)
            ic = copy(this.fast_seg);
            switch tissue
                case {'csf', 0}
                    ic = ic.numeq(1);
                    ic.fileprefix = strcat(this.fast_seg.fileprefix, 'csf');
                case {'gray', 1}
                    ic = ic.numeq(2);
                    ic.fileprefix = strcat(this.fast_seg.fileprefix, 'gray');
                case {'white', 2}
                    ic = ic.numeq(3);
                    ic.fileprefix = strcat(this.fast_seg.fileprefix, 'white');
                case {'all', 'icv'}
                    ic = ic.numgt(0);
                    ic.fileprefix = strcat(this.fast_seg.fileprefix, 'icv');
                otherwise
                    error('mladni:ValueError', 'FDG.fast_seg_binary_for()')
            end
            %disp(ic)
        end
        function        finalize(this)
            deleteExisting(strcat(this.t1w.fqfp, sprintf('_b%i.*', 10*this.blur)));
            deleteExisting(strcat(this.t1w_n4.fqfp, sprintf('_b%i.*', 10*this.blur)));
            
            if ~this.debug
                deleteExisting(this.t1w);
                deleteExisting(this.fdg);
            end
        end
        function this = flirt_fdg2t1w(this)
            %% #FLIRT_FDG2T1W

            this.flirt_fdg = mlfsl.Flirt( ...
                'in', this.fdg.fqfilename, ...
                'ref', this.t1w_blurred.fqfilename, ...
                'out', this.niigz(this.fdg_on_t1w), ...
                'omat', this.mat(this.fdg_on_t1w), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'searchrx', 180, ...
                'interp', 'trilinear');
            this.flirt_fdg.flirt();

            j0 = fileread(strcat(this.fdg.fqfp, ".json"));
            [~,j1] = this.flirt_fdg.cost_final();
            jsonrecode(j0, j1, 'filenameNew', this.json(this.fdg_on_t1w));
        end
        function this = flirt_fdg2t1w2atl(this)
            %% #FLIRT_FDG2T1W2ATL

            assert(~isempty(this.flirt_t1w))
            assert(isfile(this.niigz(this.t1w_on_atl)))

            this.flirt_fdg = mlfsl.Flirt( ...
                'in', this.fdg.fqfilename, ...
                'ref', this.t1w_blurred.fqfilename, ...
                'out', this.niigz(this.fdg_on_t1w), ...
                'omat', this.mat(this.fdg_on_t1w), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'searchrx', 180, ...
                'interp', 'trilinear');
            this.flirt_fdg.flirt();

            j0 = fileread(strcat(this.fdg.fqfp, ".json"));
            [~,j1] = this.flirt_fdg.cost_final();
            jsonrecode(j0, j1, 'filenameNew', this.json(this.fdg_on_t1w));

            flirt_on_atl_ = copy(this.flirt_fdg);
            flirt_on_atl_.concatXfm('BtoC', this.flirt_t1w.omat, 'AtoC', this.mat(this.fdg_on_atl))
            flirt_on_atl_.ref = this.atl;
            flirt_on_atl_.out = this.niigz(this.fdg_on_atl);
            flirt_on_atl_.applyXfm()

            j0 = fileread(strcat(this.fdg.fqfp, ".json"));
            [~,j1] = this.flirt_fdg.cost_final();
            jsonrecode(j0, j1, 'filenameNew', this.json(this.fdg_on_atl));
        end
        function this = flirt_t1w2atl(this)
            %% #FLIRT_T1W2ATL

            % ========== DEBUGGING chpc3 ==========
            %this.t1w
            %this.atl
            %this.t1w_on_atl
            %this.atl_mask
            %this.t1w_dlicv            
            
            this.flirt_t1w = mlfsl.Flirt( ...
                'in', this.t1w.fqfilename, ...
                'ref', this.atl.fqfilename, ...
                'out', this.niigz(this.t1w_on_atl), ...
                'omat', this.mat(this.t1w_on_atl), ...
                'bins', 256, ...
                'cost', 'corratio', ...
                'dof', 12, ...
                'interp', 'trilinear', ...
                'refweight', this.niigz(this.atl_mask), ...
                'inweight', this.niigz(this.t1w_dlicv));
            this.flirt_t1w.flirt()

            j0 = fileread(strcat(this.t1w.fqfp, ".json"));
            [~,j1] = this.flirt_t1w.cost_final();
            jsonrecode(j0, j1, 'filenameNew', this.json(this.t1w_on_atl));
        end
        function icv  = icv_from_imaging(~, ic)
            %  Returns:
            %      icv:  scalar intracranial volume in L

            ic_ = copy(ic);
            ic_ = ic_.binarized();
            ifc_ = ic_.imagingFormat;
            icv = dipsum(ic_) .* (prod(ifc_.mmppix)/1e6);
        end
        function fqfn = in_derivatives_path(this, fqfn)
            [pth,fp,x] = myfileparts(fqfn);
            pth = strrep(pth, this.rawdata_path_, this.derivatives_path_);
            fqfn = fullfile(pth, fp, x);
        end
        function this = N4BiasFieldCorrection(this)
            %if contains(this.t1w_on_atl.fileprefix, '-n3') || contains(this.t1w_on_atl.fileprefix, '-n4')
            %    this.t1w_n4_ = copy(this.t1w);
            %    return
            %end
            this.ants.N4BiasFieldCorrection(this.t1w, this.t1w_n4);
            
            mladni.FDG.jsonrecode( ...
                this.t1w, ...
                struct('image_mass', this.image_mass(this.t1w_n4)), ...
                this.t1w_n4);
        end
        function this = resolve_fdg2t1w(this)
            pwd0 = pushd(this.fdg.filepath);
            
            msks{1} = this.t1w_dlicv;
            msks{2} = mlfourd.ImagingContext2('none.nii.gz');
            imgs{1} = this.t1w_blurred;
            imgs{2} = this.fdg;
            t4rb = mlfourd.SimpleT4ResolveBuilder( ...
                'workpath', this.fdg.filepath, ...
                'maskForImages', msks, ...
                'theImages', imgs, ...
                'debug', this.debug);
            t4rb = t4rb.resolve();
            
            movefile(this.niigz(t4rb.theImagesFinal{2}), this.niigz(this.fdg_on_t1w));
            movefile(this.json(t4rb.theImagesFinal{2}), this.json(this.fdg_on_t1w));    
            mladni.FDG.jsonrecode( ...
                this.fdg_on_t1w, ...
                struct('image_activity', this.image_mass(this.fdg_on_t1w)), ...
                this.fdg_on_t1w);
            
            if ~this.debug
                t4rb.deleteFourdfp(t4rb.theImages);
                t4rb.deleteFourdfp(t4rb.theImagesOp(:,1));
            end

            popd(pwd0);
        end
        function this = save_qc(this)
            import mladni.FDGQC.*;
            
            % DEBUG
            %disp(this.atl)
            %disp(this.t1w_warped)
            %disp(this.fast_gray_detJ_warped)
            %disp(this.fdg_detJ_warped)
            
            try
                this.atl.save_qc(this.t1w_warped);
                this.t1w_warped.save_qc(this.fast_csf_detJ_warped)
                this.t1w_warped.save_qc(this.fast_gray_detJ_warped)
                this.t1w_warped.save_qc(this.fast_white_detJ_warped)
                this.t1w_warped.save_qc(this.fdg_detJ_warped)
            catch
            end
        end
        function s    = suvr_pons_vermis(this)
            t = this.table_fdg1();
            
            %disp(this.Subject(this.fdg_on_t1w))
            %disp(this.AcqDate(this.fdg_on_t1w))
            
            select = contains(t.Subject, this.Subject(this.fdg_on_t1w));
            select1 = contains(t.Description, 'Coreg, Avg, Standardized Image and Voxel Size');
            select2 = t.AcqDate == this.AcqDate(this.fdg_on_t1w);

            s = t{select & select1 & select2, 'PonsVermis'};
            s = s(~isnan(s));
            s = s(1);
        end
        function t    = table_fdg1(this)
            if ~isempty(this.table_fdg1_)
                t = this.table_fdg1_;
                return
            end
            this.table_fdg1_ = this.demog_.adni_demographics.table_fdg1();
            t = this.table_fdg1_;
        end
        function [s,r] = view(this)
            cmd = sprintf('fsleyes %s %s %s %s', this.fdg_detJ_warped.fqfn, this.fast_gray_detJ_warped.fqfn, this.t1w_warped.fqfn, this.atl.fqfn);
            [s,r] = mlbash(cmd);            
        end

        function this = FDG(varargin)
            %% #FDG 
            %  Args:
            %      fdg (required any): understood by mlfourd.ImagingContext2.  rawdata gets copied to derivatives.
            %      t1w (any): understood by mlfourd.ImagingContext2.
            %      atl (any): param understood by mlfourd.ImagingContext2.
            %      rawdata_path (text): path to replace path containing sub-*/ses-*/{pet,anat}.  Default := fdg.
            %      derivatives_path (text): path to replace path containing sub-*/ses-*/pet.  Default := fdg.
            %      blur (scalar): default := 0.
            %      proc (text): default := 'CASU'.
            
            ip = inputParser;
            addRequired(ip, "fdg")
            addParameter(ip, "t1w", [])
            addParameter(ip, "atl", fullfile(getenv("REFDIR"), "MNI152_T1_2mm.nii.gz"))
            addParameter(ip, "rawdata_path", "", @istext)
            addParameter(ip, "derivatives_path", "", @istext)
            addParameter(ip, "blur", 0, @isscalar)
            addParameter(ip, "proc", "CASU", @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;   
            ic__ = mlfourd.ImagingContext2(ipr.fdg);
            deriv_pth = strrep(ic__.filepath, 'rawdata', 'derivatives');         
            
            % construct fdg
            this.fdg_proc_ = ipr.proc;
            this.fdg_ = this.prepare_derivatives(mlfourd.ImagingContext2(ipr.fdg), deriv_pth);

            % construct blur
            this.blur_ = ipr.blur;
            if contains(this.fdg_.fileprefix, 'CASU') 
                % match blur; https://adni.loni.usc.edu/methods/pet-analysis-method/pet-analysis/
                this.blur_ = 8; % mm FWHM
            end
            
            % construct t1w
            if isempty(ipr.t1w) % find T1w nearest in time
                ipr.t1w = this.findT1w();
            end
            if isempty(ipr.t1w)
                this.t1w_ = [];
                return
            end
            this.t1w_ = this.prepare_derivatives(mlfourd.ImagingContext2(ipr.t1w), deriv_pth);

            % construct atl
            this.atl_ = mlfourd.ImagingContext2(ipr.atl);

            % construct paths
            if 0 == strlength(ipr.rawdata_path)
                re = regexp(this.fdg_.filepath, '(?<rdp>\S+)/sub-\d{3}S\d{4}/ses-\d+/pet', 'names');
                this.rawdata_path_ = re.rdp;
                assert(isfolder(this.rawdata_path_))
            end
            if 0 == strlength(ipr.derivatives_path)
                this.derivatives_path_ = fullfile(fileparts(this.rawdata_path_), 'derivatives', '');
                assert(~isempty(this.derivatives_path_))
            end

            % ANTs
            this.ants = mlfsl.ANTs('workpath', this.t1w_.filepath);
            
            % AdniDemographics
            this.demog_ = mladni.FDGDemographics();
            this.fdg_reference_ = 'ADNI';
        end
    end

    %% PROTECTED

    properties (Access = protected)
        blur_
        
        atl_
        atl_mask_
        atl_mask1_
        
        fast_csf_
        fast_csf_detJ_warped_
        fast_csf_warped_
        fast_gray_
        fast_gray_detJ_warped_        
        fast_gray_warped_
        fast_seg_
        fast_white_
        fast_white_detJ_warped_
        fast_white_warped_
        
        fdg_
        fdg_mskt_
        fdg_on_atl_
        fdg_on_t1w_
        fdg_detJ_warped_
        fdg_detJ_warped_mask_
        fdg_proc_
        fdg_reference_
        fdg_warped_
        fdg_warped_dlicv_
        
        table_fdg1_
        
        t1w_
        t1w_brain_
        t1w_blurred_
        t1w_detJ_
        t1w_dlicv_   
        t1w_dlicv_detJ_warped_     
        t1w_mask_
        t1w_n4_
        t1w_on_atl_
        t1w_warped_
        t1w_0genericaffine_
        t1w_1warp_
        t1w_2warp_

        derivatives_path_
        demog_
        rawdata_path_
    end

    methods (Access = protected)
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.fdg_ = copy(this.fdg_);
            that.t1w_ = copy(this.t1w_);
            that.atl_ = copy(this.atl_);
            that.fdg_mskt_ = copy(this.fdg_mskt_);
            that.t1w_mask_ = copy(this.t1w_mask_);
            that.atl_mask_ = copy(this.atl_mask_);
            that.fdg_on_t1w_ = copy(this.fdg_on_t1w_);
            that.t1w_on_atl_ = copy(this.t1w_on_atl_);
            that.fdg_on_atl_ = copy(this.fdg_on_atl_);
        end
        function fqfn = findT1w(this)
            fqfn = '';
            
            try
                % find /pth/to/rawdata/sub-123S4566/ses-yyyymmdd/anat/*-T1w.nii.gz
                fdgpth_ = strrep(this.fdg.filepath, 'derivatives', 'rawdata'); % /pth/to/rawdata/sub-123S4566/ses-yyyymmdd/pet
                sesfold_ = mybasename(myfileparts(fdgpth_)); % ses-yyyymmdd
                re = regexp(sesfold_, 'ses-(?<dt>\d{8})', 'names');
                dt_fdg_ = datetime(re.dt, 'InputFormat', 'yyyyMMdd');

                subpth_ = fileparts(myfileparts(fdgpth_)); % /pth/to/rawdata/sub-123S4566
                globt1w_ = globT(fullfile(subpth_, 'ses-*', 'anat', '*_T1w.nii.gz'));
                globt1w_ = globt1w_(~contains(globt1w_, 'mask'));
                if any(contains(globt1w_, 'irspgr')) % KLUDGE for irspgr containing t1w and t1w containing localizer                    
                    globt1w_ = globt1w_(contains(globt1w_, 'irspgr'));
                end
                if isempty(globt1w_)
                    return
                end                
                
                % find globt1w__ for T1w datetime closest to FDG datetime
                dts_ = NaT(size(globt1w_));
                for ig = 1:length(globt1w_)
                    sesfold_ = mybasename(myfileparts(myfileparts(globt1w_{ig})));
                    re = regexp(sesfold_, 'ses-(?<dt>\d{8})', 'names');
                    dts_(ig) = datetime(re.dt, 'InputFormat', 'yyyyMMdd');
                end
                if isempty(dts_)
                    return
                end
                ddt = min(abs(dts_ - dt_fdg_)); % ddt of anat nearest in time to fdg
                select = abs(dts_ - dt_fdg_) == ddt;
                globt1w__ = globt1w_(select);
                
                % find idx_proc of proc_ with longest processing descriptor
                proc_ = cell(size(globt1w__));
                for ip = 1:length(globt1w__)
                    [~,fp_] = myfileparts(globt1w__{ip});
                    re = regexp(fp_, 'sub-\d{3}S\d{4}_ses-\d+_acq-\w+_proc-(?<proc>\S+)_T1w', 'names');
                    proc_{ip} = re.proc;
                end
                [~,idx_proc] = max(cell2mat(cellfun(@length , proc_, 'UniformOutput', false))); % idx of longest proc specification
                if isempty(idx_proc)
                    return
                end

                % ========== DEBUGGING chpc3 ==========
                % fdgpth_
                % sesfold_
                % dt_fdg_
                % subpth_
                % fprintf('mladni.FDG.findT1w.globpth_:\n')
                % disp(ascol(globpth_))
                % dts_'
                % idx_anat
                % fprintf('mladni.FDG.findT1w.globt1w_:\n');
                % disp(ascol(globt1w_))
                % fprintf('mladni.FDG.findT1w.proc_:\n');
                % disp(ascol(proc_))                
                % idx_proc

                % find globt1w__{idx_proc}
                fqfn = globt1w__{idx_proc};
                if ~isfile(fqfn)
                    return
                end
                if ~isempty(getenv('DEBUG'))
                    fprintf('mladn.FDG.findT1w:  found %s\n', fqfn);
                end
            catch ME
                handwarning(ME)
            end
        end
        function icd  = prepare_derivatives(~, ic, deriv_pth)
            if ~isfolder(deriv_pth)
                mkdir(deriv_pth);
            end
            try
                if ~isfile(fullfile(deriv_pth, ic.filename))
                    mlbash(sprintf('cp -f %s %s', ic.fqfn, deriv_pth), 'echo', true);
                end
                if ~isfile(fullfile(deriv_pth, strcat(ic.fileprefix, '.json'))) && ...
                        isfile(strcat(ic.fqfp, '.json'))
                    mlbash(sprintf('cp -f %s %s', strcat(ic.fqfp, '.json'), deriv_pth), 'echo', true);
                end
            catch
            end
            mlbash(sprintf('chmod -R 755 %s', deriv_pth));

            icd = mlfourd.ImagingContext2(fullfile(deriv_pth, ic.filename));
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
