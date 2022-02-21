classdef FDG
    %% FDG
    %  
    %  Created 31-Dec-2021 00:02:37 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.11.0.1809720 (R2021b) Update 1 for MACI64.  Copyright 2021 John J. Lee.

    properties 
        flirt_fdg
        flirt_t1w

        fdg % mlfourd.ImagingContext2
        t1w % mlfourd.ImagingContext2
        atl % mlfourd.ImagingContext2
    end

    properties (Dependent)
        fdg_on_t1w_niigz
        fdg_on_t1w_mat
        fdg_on_t1w_json

        fdg_on_atl_niigz
        fdg_on_atl_mat
        fdg_on_atl_json

        t1w_on_atl_niigz
        t1w_on_atl_mat
        t1w_on_atl_json
    end

    methods

        %% GET
        
        function g = get.fdg_on_t1w_niigz(this)
            g = strcat(this.fdg.fqfp, "_on_T1w.nii.gz");
            g = this.in_derivatives_path(g);
        end
        function g = get.fdg_on_t1w_mat(this)
            fqfp = extractBefore(this.fdg_on_t1w_niigz, '.nii.gz');
            g = strcat(fqfp, '.mat');
            g = this.in_derivatives_path(g);
        end
        function g = get.fdg_on_t1w_json(this)
            fqfp = extractBefore(this.fdg_on_t1w_niigz, '.nii.gz');
            g = strcat(fqfp, '.json');
            g = this.in_derivatives_path(g);
        end

        function g = get.fdg_on_atl_niigz(this)
            g = strcat(this.fdg.fqfp, "_on_atl.nii.gz");
            g = this.in_derivatives_path(g);
        end
        function g = get.fdg_on_atl_mat(this)
            fqfp = extractBefore(this.fdg_on_atl_niigz, '.nii.gz');
            g = strcat(fqfp, '.mat');
            g = this.in_derivatives_path(g);
        end
        function g = get.fdg_on_atl_json(this)
            fqfp = extractBefore(this.fdg_on_atl_niigz, '.nii.gz');
            g = strcat(fqfp, '.json');
            g = this.in_derivatives_path(g);
        end

        function g = get.t1w_on_atl_niigz(this)
            g = strcat(this.t1w.fqfp, "_on_atl.nii.gz");
            g = this.in_derivatives_path(g);
        end
        function g = get.t1w_on_atl_mat(this)
            fqfp = extractBefore(this.t1w_on_atl_niigz, '.nii.gz');
            g = strcat(fqfp, '.mat');
            g = this.in_derivatives_path(g);
        end
        function g = get.t1w_on_atl_json(this)
            fqfp = extractBefore(this.t1w_on_atl_niigz, '.nii.gz');
            g = strcat(fqfp, '.json');
            g = this.in_derivatives_path(g);
        end

        %%

        function this = flirt(this)
            this = this.flirt_t1w2atl();
            this = this.flirt_fdg2t1w2atl();
        end
        function [s,r] = view(this)
            cmd = sprintf('fsleyes %s %s %s', this.fdg_on_atl_niigz, this.t1w_on_atl_niigz, this.atl.fqfn);
            [s,r] = mlbash(cmd);            
        end

        function this = flirt_fdg2t1w(this)
            this.flirt_fdg = mlfsl.Flirt( ...
                'in', this.fdg.fqfilename, ...
                'ref', this.t1w.fqfilename, ...
                'out', this.fdg_on_t1w_niigz, ...
                'omat', this.fdg_on_t1w_mat, ...
                'bins', 256, ...
                'cost', 'normmi', ...
                'dof', 6, ...
                'interp', 'trilinear');
            this.flirt_fdg.flirt()

            j0 = fileread(strcat(this.fdg.fqfp, ".json"));
            [~,j1] = this.flirt_fdg.cost_final();
            j = jsonrecode(j0, j1, 'filenameNew', this.fdg_on_t1w_json);
        end
        function this = flirt_fdg2t1w2atl(this)

            assert(~isempty(this.flirt_t1w))
            assert(isfile(this.t1w_on_atl_niigz))

            this.flirt_fdg = mlfsl.Flirt( ...
                'in', this.fdg.fqfilename, ...
                'ref', this.t1w.fqfilename, ...
                'out', this.fdg_on_t1w_niigz, ...
                'omat', this.fdg_on_t1w_mat, ...
                'bins', 256, ...
                'cost', 'normmi', ...
                'dof', 6, ...
                'interp', 'trilinear');
            this.flirt_fdg.flirt()

            j0 = fileread(strcat(this.fdg.fqfp, ".json"));
            [~,j1] = this.flirt_fdg.cost_final();
            j = jsonrecode(j0, j1, 'filenameNew', this.fdg_on_t1w_json);

            flirt_on_atl_ = copy(this.flirt_fdg);
            flirt_on_atl_.concatXfm('BtoC', this.flirt_t1w.omat, 'AtoC', this.fdg_on_atl_mat)
            flirt_on_atl_.ref = this.atl;
            flirt_on_atl_.out = this.fdg_on_atl_niigz;
            flirt_on_atl_.applyXfm

            j0 = fileread(strcat(this.fdg.fqfp, ".json"));
            [~,j1] = this.flirt_fdg.cost_final();
            j = jsonrecode(j0, j1, 'filenameNew', this.fdg_on_atl_json);
        end
        function this = flirt_t1w2atl(this)
            this.flirt_t1w = mlfsl.Flirt( ...
                'in', this.t1w.fqfilename, ...
                'ref', this.atl.fqfilename, ...
                'out', this.t1w_on_atl_niigz, ...
                'omat', this.t1w_on_atl_mat, ...
                'bins', 256, ...
                'cost', 'corratio', ...
                'dof', 12, ...
                'interp', 'trilinear');
            this.flirt_t1w.flirt()

            j0 = fileread(strcat(this.t1w.fqfp, ".json"));
            [~,j1] = this.flirt_t1w.cost_final();
            j = jsonrecode(j0, j1, 'filenameNew', this.t1w_on_atl_json);
        end
        function fqfn = in_derivatives_path(this, fqfn)
            fqfn = strrep(fqfn, this.rawdata_path_, this.derivatives_path_);
        end

        function this = FDG(varargin)
            %% FDG 
            %  Args:
            %      fdg (any): understood by mlfourd.ImagingContext2.
            %      t1w (any): understood by mlfourd.ImagingContext2.
            %      atl (any): param understood by mlfourd.ImagingContext2.
            %      derivatives_path (text): path to replace path containing sub-*/ses-*/{pet,anat}.
            
            ip = inputParser;
            addRequired(ip, "fdg")
            addRequired(ip, "t1w")
            addParameter(ip, "atl", fullfile(getenv("FSLDIR"), "data", "standard", "MNI152_T1_1mm.nii.gz"))
            addParameter(ip, "rawdata_path", "", @istext)
            addParameter(ip, "derivatives_path", "", @istext)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.fdg = mlfourd.ImagingContext2(ipr.fdg);
            this.t1w = mlfourd.ImagingContext2(ipr.t1w);
            this.atl = mlfourd.ImagingContext2(ipr.atl);

            if 0 == strlength(ipr.rawdata_path)
                re = regexp(this.fdg.filepath, '(?<rdp>\S+)/sub-\d{3}S\d{4}/ses-\d+/pet', 'names');
                this.rawdata_path_ = re.rdp;
                assert(isfolder(this.rawdata_path_))
            end
            if 0 == strlength(ipr.derivatives_path)
                this.derivatives_path_ = fullfile(fileparts(this.rawdata_path_), 'derivatives', '');
                assert(~isempty(this.derivatives_path_))
            end
        end
    end

    %% PROTECTED

    properties (Access = protected)
        derivatives_path_
        rawdata_path_
    end

    methods (Access = protected)
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
