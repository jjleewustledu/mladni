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
            %
            % >> mladni.FDG.getDebugLog(j, c)
            % LOG FILE OUTPUT:
            % --------------------------------------------------------------
            % Begin Slurm Prologue Sun Apr 24 21:44:02 CDT 2022 1650854642
            % Job ID:		1245566
            % Username:	jjlee
            % Partition:	small
            % End Slurm Prologue Sun Apr 24 21:44:02 CDT 2022 1650854642
            % --------------------------------------------------------------
            % The scheduler has allocated the following nodes to this job:
            % node15
            % "/export/matlab/R2021a//bin/mw_mpiexec" -l -n 3 "/export/matlab/R2021a/bin/worker" -parallel
            % [0] Sending a stop signal to all the labs...
            % [0] Parallel pool is shutting down.[0] 
            % 
            % ===================================================================================
            % =   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
            % =   PID 1517883 RUNNING AT node15
            % =   EXIT CODE: 1
            % =   CLEANING UP REMAINING PROCESSES
            % =   YOU CAN IGNORE THE BELOW CLEANUP MESSAGES
            % ===================================================================================
            % YOUR APPLICATION TERMINATED WITH THE EXIT STRING: Hangup (signal 1)
            % This typically refers to a problem with your application.
            % Please see the FAQ page for debugging suggestions
            % Exiting with code: 1
            % --------------------------------------------------------------
            % Begin Slurm Epilogue Sun Apr 24 21:47:32 CDT 2022 1650854852
            % Name                : Job32
            % User                : jjlee
            % Partition           : small
            % Nodes               : node15
            % Cores               : 3
            % State               : FAILED
            % Submit              : 2022-04-24T21:44:02
            % Start               : 2022-04-24T21:44:02
            % End                 : 2022-04-24T21:47:29
            % Reserved Walltime   : 01:00:00
            % Used Walltime       : 00:03:27
            % Used CPU Time       : 00:03:37
            % % User (Computation): 94.54%
            % % System (I/O)      :  5.46%
            % Mem Reserved        : 24000M
            % Max Mem Used        : 4.32G (4635340800.0)
            % Max Disk Write      : 1.18G (1268105871.36)
            % Max Disk Read       : 3.31G (3554253209.6)
            % Max-Mem-Used Node   : node15
            % Max-Disk-Write Node : node15
            % Max-Disk-Read Node  : node15
            % End Slurm Epilogue Sun Apr 24 21:47:32 CDT 2022 1650854852
            % --------------------------------------------------------------
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.FDG.batch_tiny, 1, {'len', [3]}, 'Pool', 3, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false); %#ok<NBRAK> % {'len', []}, 'Pool', 31            
        end
        function [j,c] = parcluster()
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
            for ji = 1:length(ld.globbed)
                j{ji} = c.batch(@mladni.FDG.batch, 1, {'globbed', ld.globbed{ji}}, 'Pool', 31, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);  %#ok<AGROW>
            end
        end
        function [j,c] = parcluster2()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            %  --------------------------------------------------------------
            %
            % >> mladni.FDG.getDebugLog(j, c)
            % LOG FILE OUTPUT:
            % --------------------------------------------------------------
            % --------------------------------------------------------------
            
            c = parcluster;
            disp(c.AdditionalProperties)
            
            globbing_file = fullfile(getenv('SINGULARITY_HOME'), 'ADNI', 'bids', 'derivatives', 'globbed.mat');
            assert(isfile(globbing_file));
            ld = load(globbing_file);
            g = [ld.globbed{:}];
            g1 = cellfun(@(x) strrep(strrep(x, 'rawdata', 'derivatives'), '.nii.gz', '_on_T1w_Warped.nii.gz'), ...
                g, 'UniformOutput', false);
            select = cellfun(@(x) ~isfile(x), g1);
            j = c.batch(@mladni.FDG.batch, 1, {'globbed', g(select)}, 'Pool', 31, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false); % {'len', []}, 'Pool', 31            
        end
        function t = batch_tiny(varargin)
            %% Requires completion of mladni.AdniBidsT1w.batch, mladniAdniBidsFdg.batch,
            %  and mladni.AdniBids.adni_3dresample to populate ADNI_FDG/bids/rawdata.
            %  Params:
            %      len (numeric):  # FDG sessions
            %      proc (text):  default 'CASU'
            
            ip = inputParser;
            addParameter(ip, 'len', [1 3], @isnumeric)
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            if isscalar(ipr.len)
                ipr.len = 1:ipr.len;
            end
            
            disp('Start mladni.FDG.batch()')        
            
            globbed = globT( ...
                fullfile( ...
                    '/home/aris_data/ADNI_FDG/bids/rawdata', ...
                    'sub-*', 'ses-*', 'pet', ...
                    sprintf('*-%s_*_pet.nii.gz', ipr.proc))); % *_proc-CASU_orient-rpi_pet
            globbed = globbed(ipr.len);
            t = cell2table(globbed', 'VariableNames', {'rawdata_pet_filename'});
            writetable(t, '/home/aris_data/ADNI_FDG/bids/derivatives/globbed_tiny.csv');
            fprintf('size(mladni.FDG.batch.globbed) -> %s\n', mat2str(size(globbed)));            
            
            t0 = tic;
            mladni.CHPC3.clean_tempdir();
            parfor idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});
                    %derpth = strrep(fdg_.filepath, 'rawdata', 'derivatives');
                    %if ~isfolder(derpth) % /path/to/bids/derivatives/sub-*/ses-*/pet/
                        
                        %c{idx} = derpth; %% DEBUG
                        %continue
                        
                        this = mladni.FDG(fdg_);
                        if ~isempty(this.t1w)
                            this.call_resolve(); 
                            this.finalize(); 
                        end
                    %end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            t = toc(t0);

            disp('mladni.FDG.batch() completed')            
        end
        function globbed = batch_globbed(varargin)
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
                    %if ~isfolder(derpth) % /path/to/bids/derivatives/sub-*/ses-*/pet/

                    %c{idx} = derpth; %% DEBUG
                    %continue

                    this = mladni.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.call_resolve(); 
                        this.finalize(); 
                    end
                    %end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end            
            %save('/scratch/jjlee/Singularity/ADNI/bids/derivatives/c.mat', 'c'); %% DEBUG            
            t = toc(t0);

            disp('mladni.FDG.batch() completed')            
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
        function m = image_mass(ic)
            dV = voxelVolume(ic);
            ic1 = ic.thresh(0);
            sumDensities = dipsum(ic1);            
            m = sumDensities*dV;
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
    end

    properties
        ants
        flirt_fdg
        flirt_t1w
    end

    properties (Dependent)
        blur % for co-registration of t1w with fdg
        debug % ~isempty(getenv('DEBUG'))

        atl
        atl_brain % skull-stripped brain only
        atl_mask 
        atl_mask1 
        
        fast_csf
        fast_csf_detJ_warped
        fast_gray
        fast_gray_detJ_warped
        fast_white
        fast_white_detJ_warped
        
        fdg % mlfourd.ImagingContext2
        fdg_mask % mlfourd.ImagingContext2 
        fdg_mskt
        fdg_on_atl
        fdg_on_t1w % mlfourd.ImagingContext2 
        fdg_detJ_warped % warped to atl, weighted by det(J)
        fdg_detJ_warped_mask % generated from dlicv, binarized, 8mm blurring, thresh 0.1, binarized
        
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

        %% GET

        function g = get.blur(this)
            g = this.blur_;
        end
        function g = get.debug(~)
            g = ~isempty(getenv('DEBUG'));
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
        function this = build_fast_warped(this)
            % warp all segs to atlas    
            % Args:
            %     in (any): pve from fast understood by ImagingContext2.
            
            for pve = {'fast_csf', 'fast_gray', 'fast_white'}                
                this.([pve{1} '_detJ_warped_']) = this.build_fast_warped_seg(this.(pve{1}));
            end
        end
        function outw = build_fast_warped_seg(this, in)
            inw = mlfourd.ImagingContext2( ...
                this.ants.antsApplyTransforms(in, this.atl_brain, this.t1w_brain));
            outw = this.t1w_detJ .* inw;
            outw.fqfp = strcat(in.fqfp, '_detJ_Warped');
            outw.save();

            mladni.FDG.jsonrecode( ...
                inw, ...
                struct('state_changes', 'outw = this.t1w_detJ .* inw', ...
                       'image_mass', this.image_mass(outw)), ...
                outw);

            deleteExisting(inw);
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
        function this = call_resolve(this)
            %% #CALL_RESOLVE
            
            this = this.N4BiasFieldCorrection();
            this = this.deepmrseg_apply();
            this = this.resolve_fdg2t1w();
            this = this.fast();
            this = this.CreateJacobianDeterminantImages();
            this = this.build_deepmrseg_warped();
            this = this.build_fast_warped();
            this = this.build_fdg_warped();
            this = this.build_nmf_inputs();
            %this = this.save_qc();
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
            %      blur (scalar):  Default := 0.
            
            ip = inputParser;
            addRequired(ip, "fdg")
            addParameter(ip, "t1w", [])
            addParameter(ip, "atl", fullfile(getenv("REFDIR"), "MNI152_T1_2mm.nii.gz"))
            addParameter(ip, "rawdata_path", "", @istext)
            addParameter(ip, "derivatives_path", "", @istext)
            addParameter(ip, "blur", 0, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;   
            ic__ = mlfourd.ImagingContext2(ipr.fdg);
            deriv_pth = strrep(ic__.filepath, 'rawdata', 'derivatives');         
            
            % construct fdg
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
        fast_gray_
        fast_gray_detJ_warped_
        fast_white_
        fast_white_detJ_warped_
        
        fdg_
        fdg_mskt_
        fdg_on_atl_
        fdg_on_t1w_
        fdg_detJ_warped_
        fdg_detJ_warped_mask_
        
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
                fdgpth_ = strrep(this.fdg.filepath, 'derivatives', 'rawdata'); % /pth/to/rawdata/sub-123S4566/ses-yyyymmdd/pet
                sesfold_ = mybasename(myfileparts(fdgpth_)); % ses-yyyymmdd
                re = regexp(sesfold_, 'ses-(?<dt>\d{8})', 'names');
                dt_fdg_ = datetime(re.dt, 'InputFormat', 'yyyyMMdd');

                subpth_ = fileparts(myfileparts(fdgpth_)); % /pth/to/rawdata/sub-123S4566
                globpth_ = globT(fullfile(subpth_, 'ses-*', 'anat', ''));
                if isempty(globpth_)
                    return
                end
                
                
                
                dts_ = NaT(size(globpth_));
                for ig = 1:length(globpth_)
                    sesfold_ = mybasename(myfileparts(globpth_{ig}));
                    re = regexp(sesfold_, 'ses-(?<dt>\d{8})', 'names');
                    dts_(ig) = datetime(re.dt, 'InputFormat', 'yyyyMMdd');
                end
                if isempty(dts_)
                    return
                end
                [~,idx_anat] = min(abs(dts_ - dt_fdg_)); % idx of anat nearest in time to fdg
                if isempty(idx_anat)
                    return
                end
                
                

                globt1w_ = globT(fullfile(globpth_{idx_anat}, '*_T1w.nii.gz')); % /pth/to/rawdata/sub-123S4566/ses-yyyyMMdd/anat/sub-*_ses-*_acq-*_proc-*_T1w.nii.gz
                globt1w_ = globt1w_(~contains(globt1w_, 'mask'));
                if isempty(globt1w_)
                    return
                end
                
                
                
                proc_ = cell(size(globt1w_));
                for ip = 1:length(globt1w_)
                    [~,fp_] = myfileparts(globt1w_{ip});
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

                fqfn = globt1w_{idx_proc};
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
                mlbash(sprintf('cp -f %s %s', ic.fqfn, deriv_pth), 'echo', true);
                mlbash(sprintf('cp -f %s %s', strcat(ic.fqfp, '.json'), deriv_pth), 'echo', true);
            catch
            end
            mlbash(sprintf('chmod -R 755 %s', deriv_pth));

            icd = mlfourd.ImagingContext2(fullfile(deriv_pth, ic.filename));
%             if ~contains(icd.fileprefix, '-std')
%                 tmp = icd.fqfn;
%                 icd.reorient2std();
%                 icd.selectNiftiTool();
%                 icd.save();
%                 deleteExisting(tmp);
%             end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
