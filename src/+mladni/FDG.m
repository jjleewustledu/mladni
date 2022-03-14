classdef FDG < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% FDG uses object-oriented mlfsl.FSL for registration and mlfourd.ImagingContext2 for actions on imaging data.  
    %  Most common use case:
    %  >> an_fdg = mlfourd.ImagingContext2(fq_filename)
    %  >> obj = mladni.FDG(an_fdg);
    %  >> obj.call()
    %  >> obj.view()
    %  
    %  Created 31-Dec-2021 00:02:37 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1809720 (R2021b) Update 1 for MACI64.  Copyright 2021 John J. Lee.

    methods (Static)
        function fn = niigz(obj)
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.nii.gz');
        end
        function fn = mat(obj)
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.mat');
        end
        function fn = json(obj)
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.json');
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

        fdg % mlfourd.ImagingContext2        
        t1w  
        t1wb
        atl

        fdg_mask % mlfourd.ImagingContext2 
        t1w_mask
        atl_mask 
        atl_mask1 

        fdg_on_t1w % mlfourd.ImagingContext2 
        t1w_on_atl
        t1w_on_atl_n4
        t1w_on_atl_dmrs
        t1w_on_atl_warped
        t1w_on_atl_1warp
        t1w_on_atl_detJ
        fdg_on_atl
        fdg_final % affine registered to atl, warped to atl, weighted by det(J)
    end

    methods

        %% GET

        function g = get.blur(this)
            g = this.blur_;
        end
        function g = get.debug(~)
            g = ~isempty(getenv('DEBUG'));
        end
        
        function g = get.fdg(this)
            g = this.fdg_;
        end
        function g = get.t1w(this)
            g = this.t1w_;
        end
        function g = get.t1wb(this)
            g = this.t1wb_;
            if ~isfile(g.fqfn)
                g.save();
            end
        end
        function g = get.atl(this)
            g = this.atl_;
        end

        function g = get.fdg_mask(this)
            if ~isempty(this.fdg_mask_)
                g = this.fdg_mask_;
                return
            end
            pwd0 = pushd(this.fdg.filepath);
            this.fdg_mask_ = this.fdg.thresh(0.001*dipmax(this.fdg));
            this.fdg_mask_ = this.fdg_mask_.binarized();
            this.fdg_mask_.fileprefix = strcat(this.fdg.fileprefix, '_msk');
            this.fdg_mask_.save();
            g = this.fdg_mask_;
            popd(pwd0);
        end
        function g = get.t1w_mask(this)
            if ~isempty(this.t1w_mask_)
                g = this.t1w_mask_;
                return
            end   
            pwd0 = pushd(this.t1w.filepath);         
            this.t1w_mask_ = mlfourd.ImagingContext2(mlfsl.Flirt.msktgen(this.t1w.fqfn));
            g = this.t1w_mask_;
            popd(pwd0);
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
        
        function g = get.fdg_on_t1w(this)
            if ~isempty(this.fdg_on_t1w_)
                g = this.fdg_on_t1w_;
                return
            end
            fqfp = this.in_derivatives_path(this.fdg.fqfp);
            this.fdg_on_t1w_ = mlfourd.ImagingContext2(strcat(fqfp, '_on_T1w.nii.gz'));
            g = this.fdg_on_t1w_;
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
        function g = get.t1w_on_atl_n4(this)
            if ~isempty(this.t1w_on_atl_n4_)
                g = this.t1w_on_atl_n4_;
                return
            end
            fqfp = this.t1w_on_atl.fqfp;
            this.t1w_on_atl_n4_ = mlfourd.ImagingContext2(strrpe(fqfp, '_T1w', '-n4_T1w'));
            g = this.t1w_on_atl_n4_;
        end
        function g = get.t1w_on_atl_dmrs(this)
            if ~isempty(this.t1w_on_atl_dmrs_)
                g = this.t1w_on_atl_dmrs_;
                return
            end
            fqfp = this.t1w_on_atl.fqfp;
            this.t1w_on_atl_dmrs_ = mlfourd.ImagingContext2(strcat(fqfp, '_dmrs.nii.gz'));
            g = this.t1w_on_atl_dmrs_;
        end
        function g = get.t1w_on_atl_warped(this)
            if ~isempty(this.t1w_on_atl_warped_)
                g = this.t1w_on_atl_warped_;
                return
            end
            fqfp = this.t1w_on_atl.fqfp;
            this.t1w_on_atl_warped_ = mlfourd.ImagingContext2(strcat(fqfp, '_Warped.nii.gz'));
            g = this.t1w_on_atl_warped_;
        end
        function g = get.t1w_on_atl_1warp(this)
            if ~isempty(this.t1w_on_atl_1warp_)
                g = this.t1w_on_atl_1warp_;
                return
            end
            fqfp = this.t1w_on_atl.fqfp;
            this.t1w_on_atl_1warp_ = mlfourd.ImagingContext2(strcat(fqfp, '_1Warp.nii.gz'));
            g = this.t1w_on_atl_1warp_;
        end
        function g = get.t1w_on_atl_detJ(this)
            if ~isempty(this.t1w_on_atl_detJ_)
                g = this.t1w_on_atl_detJ_;
                return
            end
            fqfp = this.t1w_on_atl.fqfp;
            this.t1w_on_atl_detJ_ = mlfourd.ImagingContext2(strcat(fqfp, '_detJ.nii.gz'));
            g = this.t1w_on_atl_detJ_;
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
        function g = get.fdg_final(this)
            if ~isempty(this.fdg_final_)
                g = this.fdg_final_;
                return
            end
            fqfp = this.in_derivatives_path(this.fdg.fqfp);
            this.fdg_final_ = mlfourd.ImagingContext2(strcat(fqfp, '_final.nii.gz'));
            g = this.fdg_final_;
        end

        %%

        function this = call(this)
            this = this.flirt_t1w2atl();
            this = this.flirt_fdg2t1w2atl();
            this = this.N4BiasFieldCorrection();
            this = this.CreateJacobianDeterminantImage();
        end
        function this = CreateJacobianDeterminantImage(this)
            this.t1w_on_atl_warped_ = mlfourd.ImagingContext2( ...
                this.ants.antsRegistrationSyNQuick(this.atl, this.t1w_on_atl_n4));
            fdg_warped = mlfourd.ImagingContext2( ...
                this.ants.antsApplyTransforms(this.fdg_on_atl, this.atl, this.t1w_on_atl_n4));
            this.ants.CreateJacobianDeterminantImage(this.t1w_on_atl_1warp, this.t1w_on_atl_detJ);
            this.fdg_final_ = this.atl_mask1 .* this.t1w_on_atl_detJ .* fdg_warped;
            this.fdg_final_.addLog(char(this.t1w_on_atl_detJ.logger));
            this.fdg_final_.addLog(char(fdg_warped));
            this.fdg_final_.fqfileprefix = strcat(this.fdg.fqfileprefix, '_final');
            this.fdg_final_.save();
            copyfile(strcat(fdg_warped.fqfp, '.json'), strcat(this.fdg_final_.fqfp, '.json'))

            deleteExisting(fdg_warped);
        end
        function this = deepmrseg_apply(this)
        end
        function finalize(this)
            deleteExisting(strcat(this.t1w.fqfp, sprintf('_b%i.*', 10*this.blur)));
            deleteExisting(strcat(this.t1w_on_atl_n4.fqfp, '_InverseWarped.*'));

            if ~this.debug
                deleteExisting(this.fdg);
                deleteExisting(this.fdg_mask);
                deleteExisting(this.fdg_on_t1w);
                deleteExisting(this.fdg_on_atl);
                deleteExisting(this.t1w);
                deleteExisting(this.t1w_mask);
                deleteExisting(this.t1w_on_atl);
                %deleteExisting(this.t1w_on_atl_detJ);
            end
        end
        function this = flirt_fdg2t1w(this)

            this.flirt_fdg = mlfsl.Flirt( ...
                'in', this.fdg.fqfilename, ...
                'ref', this.t1wb.fqfilename, ...
                'out', this.niigz(this.fdg_on_t1w), ...
                'omat', this.mat(this.fdg_on_t1w), ...
                'bins', 256, ...
                'cost', 'normmi', ...
                'dof', 6, ...
                'interp', 'trilinear', ...
                'refweight', this.niigz(this.t1w_mask), ...
                'inweight', this.niigz(this.fdg_mask));
            this.flirt_fdg.flirt()

            j0 = fileread(strcat(this.fdg.fqfp, ".json"));
            [~,j1] = this.flirt_fdg.cost_final();
            jsonrecode(j0, j1, 'filenameNew', this.json(this.fdg_on_t1w));
        end
        function this = flirt_fdg2t1w2atl(this)

            assert(~isempty(this.flirt_t1w))
            assert(isfile(this.niigz(this.t1w_on_atl)))

            this.flirt_fdg = mlfsl.Flirt( ...
                'in', this.fdg.fqfilename, ...
                'ref', this.t1wb.fqfilename, ...
                'out', this.niigz(this.fdg_on_t1w), ...
                'omat', this.mat(this.fdg_on_t1w), ...
                'bins', 256, ...
                'cost', 'normmi', ...
                'dof', 6, ...
                'interp', 'trilinear', ...
                'refweight', this.niigz(this.t1w_mask), ...
                'inweight', this.niigz(this.fdg_mask));
            this.flirt_fdg.flirt()

            j0 = fileread(strcat(this.fdg.fqfp, ".json"));
            [~,j1] = this.flirt_fdg.cost_final();
            jsonrecode(j0, j1, 'filenameNew', this.json(this.fdg_on_t1w));

            flirt_on_atl_ = copy(this.flirt_fdg);
            flirt_on_atl_.concatXfm('BtoC', this.flirt_t1w.omat, 'AtoC', this.mat(this.fdg_on_atl))
            flirt_on_atl_.ref = this.atl;
            flirt_on_atl_.out = this.niigz(this.fdg_on_atl);
            flirt_on_atl_.applyXfm

            j0 = fileread(strcat(this.fdg.fqfp, ".json"));
            [~,j1] = this.flirt_fdg.cost_final();
            jsonrecode(j0, j1, 'filenameNew', this.json(this.fdg_on_atl));
        end
        function this = flirt_t1w2atl(this)

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
                'inweight', this.niigz(this.t1w_mask));
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
            if contains(this.t1w_on_atl.fileprefix, '-n3') || contains(this.t1w_on_atl.fileprefix, '-n4')
                this.t1w_on_atl_n4_ = copy(this.t1w_on_atl);
                return
            end
            this.ants.N4BiasFieldCorrection(this.t1w_on_atl, this.t1w_on_atl_n4);
        end        
        function [s,r] = view(this)
            cmd = sprintf('fsleyes %s %s %s', this.fdg_final.fqfn, this.t1w_on_atl_warped.fqfn, this.atl.fqfn);
            [s,r] = mlbash(cmd);            
        end

        function this = FDG(varargin)
            %% FDG 
            %  Args:
            %      fdg (required any): understood by mlfourd.ImagingContext2.  rawdata gets copied to derivatives.
            %      t1w (any): understood by mlfourd.ImagingContext2.
            %      atl (any): param understood by mlfourd.ImagingContext2.
            %      rawdata_path (text): path to replace path containing sub-*/ses-*/{pet,anat}.  Default := fdg.
            %      derivatives_path (text): path to replace path containing sub-*/ses-*/{pet,anat}.  Default := fdg.
            
            ip = inputParser;
            addRequired(ip, "fdg")
            addParameter(ip, "t1w", [])
            addParameter(ip, "atl", fullfile(getenv("FSLDIR"), "data", "standard", "MNI152_T1_2mm.nii.gz"))
            addParameter(ip, "rawdata_path", "", @istext)
            addParameter(ip, "derivatives_path", "", @istext)
            addParameter(ip, "blur", 3, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            % construct fdg
            this.fdg_ = this.prepare_derivatives(mlfourd.ImagingContext2(ipr.fdg));

            % construct t1w
            if isempty(ipr.t1w) % find T1w nearest in time
                ipr.t1w = this.findT1w();
            end
            this.t1w_ = this.prepare_derivatives(mlfourd.ImagingContext2(ipr.t1w));
            this.t1wb_ = this.t1w_.blurred(this.blur);

            % construct atl
            this.atl_ = mlfourd.ImagingContext2(ipr.atl);

            % construct blur
            this.blur_ = ipr.blur;
            if contains(this.fdg_.fileprefix, 'CASU') % match blur
                this.blur_ = 7;
            end

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
        fdg_
        t1w_
        t1wb_
        atl_
        
        fdg_mask_
        t1w_mask_
        atl_mask_
        atl_mask1_

        fdg_on_t1w_
        t1w_on_atl_
        t1w_on_atl_n4_
        t1w_on_atl_dmrs_
        t1w_on_atl_warped_
        t1w_on_atl_1warp_
        t1w_on_atl_detJ_
        fdg_on_atl_
        fdg_final_

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
            that.fdg_mask_ = copy(this.fdg_mask_);
            that.t1w_mask_ = copy(this.t1w_mask_);
            that.atl_mask_ = copy(this.atl_mask_);
            that.fdg_on_t1w_ = copy(this.fdg_on_t1w_);
            that.t1w_on_atl_ = copy(this.t1w_on_atl_);
            that.fdg_on_atl_ = copy(this.fdg_on_atl_);
        end
        function fqfn = findT1w(this)
            fdgpth_ = strrep(this.fdg.filepath, 'derivatives', 'rawdata'); % /pth/to/rawdata/sub-123S4566/ses-yyyymmdd/pet

            sesfold_ = mybasename(myfileparts(fdgpth_)); % ses-yyyymmdd
            re = regexp(sesfold_, 'ses-(?<dt>\d{8})', 'names');
            dt_fdg_ = datetime(re.dt, 'InputFormat', 'yyyyMMdd');

            subpth_ = fileparts(myfileparts(fdgpth_)); % /pth/to/rawdata/sub-123S4566
            globpth_ = globT(fullfile(subpth_, 'ses-*', 'anat', ''));
            dts_ = NaT(size(globpth_));
            for ig = 1:length(globpth_)
                sesfold_ = mybasename(myfileparts(globpth_{ig}));
                re = regexp(sesfold_, 'ses-(?<dt>\d{8})', 'names');
                dts_(ig) = datetime(re.dt, 'InputFormat', 'yyyyMMdd');
            end
            [~,idx_anat] = min(abs(dts_ - dt_fdg_)); % idx of anat nearest in time to fdg

            globt1w_ = globT(fullfile(globpth_{idx_anat}, '*_T1w.nii.gz')); % /pth/to/rawdata/sub-123S4566/ses-yyyyMMdd/anat/sub-*_ses-*_acq-*_proc-*_T1w.nii.gz
            globt1w_ = globt1w_(~contains(globt1w_, 'mask'));
            proc_ = cell(size(globt1w_));
            for ip = 1:length(globt1w_)
                [~,fp_] = myfileparts(globt1w_{ip});
                re = regexp(fp_, 'sub-\d{3}S\d{4}_ses-\d+_acq-\w+_proc-(?<proc>\S+)_T1w', 'names');
                proc_{ip} = re.proc;
            end
            [~,idx_proc] = max(cell2mat(cellfun(@length , proc_, 'UniformOutput', false))); % idx of longest proc specification

            fqfn = globt1w_{idx_proc};
        end
        function ic = prepare_derivatives(~, ic)
            if contains(ic.filepath, 'rawdata') % copy to derivatives
                ic.selectNiftiTool();
                ic.filepath = strrep(ic.filepath, 'rawdata', 'derivatives');
                ic.save();
            end
            ic.forceradiological();
            tmp = ic.fqfn;
            ic.reorient2std();
            deleteExisting(tmp);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
