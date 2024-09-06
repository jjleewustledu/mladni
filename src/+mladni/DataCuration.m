classdef (Abstract) DataCuration < handle
    %% line1
    %  line2
    %  
    %  Created 15-Apr-2024 22:05:34 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2537033 (R2024a) for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        demogr
        nmf
        nmfc
        nmfh
        nmfr
        N_patterns
    end

    properties (Abstract)
        data_home
        nmf_fdg_home
    end
    
    properties (Constant)
        selected_spans = mladni.NMFHierarchies.selected_spans
    end

    methods
        function this = DataCuration(opts)
            arguments
                opts.N_patterns double = mladni.NMF.N_PATTERNS
            end
            this.N_patterns = opts.N_patterns;
        end

        function build_stats_imaging(this, opts)
            %% This apex function calls create_montage, build_{mean_imaging, median_imaging}.

            arguments
                this mladni.DataCuration
                opts.inputDir {mustBeFolder} = fullfile(this.nmf_fdg_home, "baseline_cn")
                % e.g., NMF_FDG/baseline_cn, NMF_FDG/longitudinal_cdr_0.5_apos, ...
                % has nifti_files_mounted.csv visible to glob()
                opts.outputDir {mustBeTextScalar} = this.nmf_fdg_home 
            end
            assert(contains(opts.inputDir, "NMF"))
            if isemptytext(opts.outputDir)
                opts.outputDir = opts.inputDir;
            end
            
            % montage
            if contains(opts.inputDir, "baseline")
                globbed = glob(fullfile(opts.inputDir, sprintf("NumBases%i", this.N_patterns)));
                if isempty(globbed)
                    globbed = glob(fullfile(opts.inputDir, "*", sprintf("NumBases%i", this.N_patterns)));
                end
                assert(~isempty(globbed))
                this.nmf.create_montage(globbed{1});
            end

            % mean, var, std, median, iqr
            globbed = glob(fullfile(opts.inputDir, "nifti_files_mounted.csv"));
            if isempty(globbed)
                globbed = glob(fullfile(opts.inputDir, "*", "nifti_files_mounted.csv"));
            end
            assert(~isempty(globbed))
            this.build_mean_imaging(globbed{1}, outputDir=opts.outputDir);
            this.build_median_imaging(globbed{1}, outputDir=opts.outputDir);            
        end

        function build_mean_imaging(this, varargin)
            %% creates mean, var, std
            %  Args:
            %       csv {mustBeFile}
            %       outputDir {mustBeFolder}

            ip = inputParser;
            addRequired(ip, 'csv', @isfile);
            addParameter(ip, 'outputDir', pwd, @isfolder);
            parse(ip, varargin{:});
            ipr = ip.Results;            
            
            csv = readlines(ipr.csv); % string vector

            % create mean
            ic = mlfourd.ImagingContext2(char(csv(1)));
            ic.selectMatlabTool();
            fp = this.fileprefix_for_stat(ic.fileprefix, stat='mean');
            for i = 2:length(csv)
                try
                    ic = ic + mlfourd.ImagingContext2(char(csv(i)));
                    ic.fileprefix = fp;
                catch
                end
            end
            ic = ic/length(csv);
            ic.fileprefix = fp;
            ic.filepath = ipr.outputDir;
            ic.save();
            
            % create variance
            ic1 = mlfourd.ImagingContext2(char(csv(1)));
            ic1.selectMatlabTool();
            ic1 = (ic1 - ic).^2;
            fp1 = this.fileprefix_for_stat(ic.fileprefix, stat='variance', from_bids=false);
            for i = 2:length(csv)
                try
                    ic1 = ic1 + (mlfourd.ImagingContext2(char(csv(i))) - ic).^2;
                    ic1.fileprefix = fp1;
                catch
                end
            end
            ic1 = ic1/(length(csv) - 1);
            ic1.fileprefix = fp1;
            ic1.filepath = ipr.outputDir;
            ic1.save();

            % create std
            ic2 = copy(ic1);
            ic2 = sqrt(ic2);
            ic2.fileprefix = this.fileprefix_for_stat(ic1.fileprefix, stat='std', from_bids=false);
            ic2.save();
        end
        
        function fp = fileprefix_for_stat(this, fileprefix, opts)
            arguments
                this mladni.DataCuration
                fileprefix {mustBeTextScalar}
                opts.stat {mustBeTextScalar} = 'mean'
                opts.from_bids logical = true
                opts.starting_with {mustBeTextScalar} = 'all_'
            end

            if opts.from_bids
                re = regexp(fileprefix, '(?<sub>sub-\w+_)(?<ses>ses-(d|)\d+_)(?<descrip>\S+)', 'names');
                fp = sprintf('all_%s_%s', re.descrip, opts.stat);
                return
            end
            if startsWith(fileprefix, opts.starting_with)
                re = regexp(fileprefix, [opts.starting_with, '(?<descrip>\S+)_', opts.stat], 'names');
                fp = sprintf('all_%s_%s', re.descrip, opts.stat);
                return
            end
        end

        function build_median_imaging(this, varargin)
            %% creates median, iqr
            %  Args:
            %       csv {mustBeFile}
            %       outputDir {mustBeFolder}

            ip = inputParser;
            addRequired(ip, 'csv', @isfile);
            addParameter(ip, 'outputDir', pwd, @isfolder);
            parse(ip, varargin{:});
            ipr = ip.Results;                        

            csv = readlines(ipr.csv); % string vector

            % create median
            ifc = mlfourd.ImagingFormatContext2(char(csv(1)));
            fp = this.fileprefix_for_stat(ifc.fileprefix, stat='median'); 
            size3d = size(ifc.img);
            for i = 2:length(csv)
                try
                    ifc1 = mlfourd.ImagingFormatContext2(char(csv(i)));
                    if isempty(ifc1.img)
                        ifc.img(:,:,:,i) = nan(size3d);
                    else
                        ifc.img(:,:,:,i) = ifc1.img;
                    end
                catch
                    ifc.img(:,:,:,i) = nan(size3d);
                end
            end
            ifc.img = median(ifc.img, 4, 'omitnan');
            ifc.fileprefix = fp;
            ifc.filepath = ipr.outputDir;
            ifc.save();     

            % create iqr
            ifc = mlfourd.ImagingFormatContext2(char(csv(1)));
            fp = this.fileprefix_for_stat(ifc.fileprefix, stat='iqr');
            size3d = size(ifc.img);
            for i = 2:length(csv)
                try
                    ifc1 = mlfourd.ImagingFormatContext2(char(csv(i)));
                    if isempty(ifc1.img)
                        ifc.img(:,:,:,i) = nan(size3d);
                    else
                        ifc.img(:,:,:,i) = ifc1.img;
                    end
                catch
                    ifc.img(:,:,:,i) = nan(size3d);
                end
            end
            ifc.img = iqr(ifc.img, 4); % safe for nans
            ifc.fileprefix = fp;
            ifc.filepath = ipr.outputDir;
            ifc.save();   
        end        
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
