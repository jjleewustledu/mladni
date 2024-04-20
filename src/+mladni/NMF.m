 classdef NMF < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% baseline_cn ~ 22 bases
    %  
    %  Created 23-Jun-2022 13:02:04 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
    
    properties
        cache_files
        data_home
        downSample
        isList
        mask
        mcrroot
        memInGB
        nmfDataset
        numBases
        permute
        repetitions
        selectedNumBases
        smooth
        study_design
        volbin
    end

    properties (Constant)
        MAX_NUM_BASES = 40
        N_PATTERNS = 24
        fig_position = [80 80 1618 1000]; % coordinates for figures ~ [x0, y0, Dx, Dy]
        groups = { ...
            'cn', 'preclinical', ...
            'cdr_0p5_apos', ...
            'cdr_gt_0p5_apos', ...
            'cdr_gt_0_aneg'}
    end

    properties (Dependent)
        componentDir
        dataDemographics
        inFiles
        inFiles2
        isCrossSectional
        nmf_fdg_home
        outputDir
        targetDatasetDir % baseline_cn that provides all NMF patterns
        Xmat
        X2mat
    end

    methods % GET
        function g = get.dataDemographics(this)
            if ~isempty(this.dataDemographics_)
                g = this.dataDemographics_;
                return
            end
            if endsWith(this.data_home, 'ADNI')
                this.dataDemographics_ = mladni.AdniDemographics(study_design=this.study_design);
                g = this.dataDemographics_;
                return
            end
            if endsWith(this.data_home, 'OASIS3')
                this.dataDemographics_ = mloasis.OasisDemographics(study_design=this.study_design);
                g = this.dataDemographics_;
                return
            end
            error('mladni:RuntimeError', stackstr())
        end
        function g = get.componentDir(this)
            g = fullfile(this.outputDir, sprintf('NumBases%i', this.selectedNumBases), 'components');
            ensuredir(g);
        end
        function g = get.inFiles(this)
            g = fullfile(this.outputDir, 'nifti_files_mounted.csv');
            if ~isfile(g)
                g0 = strrep(g, '_mounted.csv', '.csv');
                t = readtable(g0, Format='%s', Delimiter=' ', ReadVariableNames=false);
                t.Var1 = strrep(t.Var1, '/scratch/jjlee', '/home/usr/jjlee');
                writetable(t, g, WriteVariableNames=false);
            end
        end
        function g = get.inFiles2(this)
            if contains(this.data_home, 'OASIS3')
                nmfc = mloasis.NMFCovariates();
                g = nmfc.inFiles;
                return
            end
            if contains(this.data_home, 'MAYO')
                nmfc = mljones.NMFCovariates();
                g = nmfc.inFiles;
                return
            end
            nmfc = mladni.NMFCovariates();
            g = nmfc.inFiles;
        end
        function g = get.isCrossSectional(this)
            g = contains(this.study_design, 'cross', IgnoreCase=true) && ...
                contains(this.study_design, 'sectional', IgnoreCase=true);
        end
        function g = get.nmf_fdg_home(this)
            g = fullfile(this.data_home, 'NMF_FDG');
        end
        function g = get.outputDir(this)
            g = fullfile(this.data_home, 'NMF_FDG', this.nmfDataset, '');
        end
        function g = get.targetDatasetDir(this)
            g = fullfile(this.data_home, 'NMF_FDG', 'baseline_cn');
        end
        function g = get.Xmat(this)
            numBasesFolder = sprintf('NumBases%i', this.selectedNumBases);
            g = fullfile(this.outputDir, numBasesFolder, 'X.mat');
        end
        function g = get.X2mat(this)
            numBasesFolder = sprintf('NumBases%i', this.selectedNumBases);
            g = fullfile(this.outputDir, numBasesFolder, 'X2.mat');  
        end
    end
    
    methods        
        function this = NMF(varargin)
            %% NMF works best with pwd ~ this.data_home/NMF_FDG
            %  Args:
            %      downSample (double): >= 1
            %      memInGB (scalar):  slurm mem request per batch job.
            %      nmfDataset (text): tag for dataset, used to specify ADNI/NMF_FDG/nmfDataset/nifti_files.csv; 
            %                         default is "test".
            %      volbin (text): home for NMF standalone executable.
            
            ip = inputParser;
            addParameter(ip, "cache_files", false, @islogical);
            addParameter(ip, "downSample", 1, @(x) x >= 1);
            addParameter(ip, "isList", true, @islogical);
            addParameter(ip, "mask", fullfile(getenv("ADNI_HOME"), "VolBin", "mask.nii.gz"), @isfile);
            addParameter(ip, "memInGB", 4, @isscalar);
            addParameter(ip, "nmfDataset", "baseline_cn", @(x) istext(x));
            addParameter(ip, "numBases", 2:2:mladni.NMF.MAX_NUM_BASES, @isnumeric)
            addParameter(ip, "permute", false, @islogical);
            addParameter(ip, "repetitions", 50, @isscalar);
            addParameter(ip, "selectedNumBases", mladni.NMF.N_PATTERNS, @isscalar);  % [2, 8, 10, 12, 14 24]
            addParameter(ip, "smooth", false, @islogical);
            addParameter(ip, "study_design", "longitudinal", @istext)
            addParameter(ip, "data_home", getenv("ADNI_HOME"), @isfolder);            
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.cache_files = ipr.cache_files;
            this.downSample = ipr.downSample;
            this.isList = ipr.isList;
            this.mask = ipr.mask;
            this.memInGB = ipr.memInGB;
            this.nmfDataset = ipr.nmfDataset;
            this.study_design = ipr.study_design;
            this.permute = ipr.permute;
            this.repetitions = ipr.repetitions;
            this.selectedNumBases = ipr.selectedNumBases;
            this.smooth = ipr.smooth;
            this.data_home = ipr.data_home;
            this.volbin = fullfile(this.data_home, "VolBin");

            assert(isfolder(this.outputDir));
        end

        function build_for_brainsmash(this)
            %% Prepares intermediates needed for brainsmash inference described in 
            %  test_neurodegeneration2.py/TestNeurodegeneration2.test_build_stats27.
            %  Requires availability of AFNI.3dmaskdump.

            arguments
                this mladni.NMF
            end            

            niiImg_home = fullfile(this.data_home, 'NMF_FDG', 'baseline_cn', ...
                sprintf('NumBases%i', this.selectedNumBases), 'OPNMF', 'niiImg');
            pwd1 = pushd(niiImg_home);
            parfor idx = 1:this.selectedNumBases
                system(sprintf( ...
                    '3dmaskdump -mask %s -o Basis_%i.txt -noijk Basis_%i.nii', ...
                    fullfile(this.data_home, 'VolBin', 'mask.nii.gz'), ...
                    idx, idx));
            end
            popd(pwd1)
        end        
        
        function rec_errors = build_repeated_calc_rec_error(this, opts)
            %% writes imaging data to X.mat, 
            %  then repeatedly calls this.calcRecError() for anticlustered bootstraps,
            %  then prepares plotting rec_errors.

            arguments
                this mladni.NMF
                opts.cache_files logical = true
                opts.subgroup {mustBeTextScalar} = 'cn'
                opts.Ncores {mustBeNumeric} = 16
            end

            param.isList = 1 ;
            param.downSample = this.downSample;
            param.smooth = 0;
            param.mask = fullfile(this.data_home, 'VolBin/mask.nii.gz');
            param.permute = 0;
            param.numBase = this.selectedNumBases;
            param.componentDir = this.componentDir;

            % BUG FIX:  selectVoxels must be logical, and
            %           numel(selectVoxels) == numel(mask_for_imaging), and
            %           size(selectVoxels) == [N_voxels 1], and
            %           size(X,1) == sum(selectVoxels)
            assert(isfile(param.mask))
            mask_ic = mlfourd.ImagingContext2(param.mask);
            selectVoxels = logical(mask_ic);
            selectVoxels = selectVoxels(:);

            resultsDir_ = fullfile(this.nmf_fdg_home, sprintf('baseline_%s', opts.subgroup), 'results');
            ensuredir(resultsDir_);

            %% build the estimates for each anticlustered bootstrap

            if ~isfile(fullfile(resultsDir_, [stackstr() '.mat']))
                pathDirs = glob( ...
                    fullfile(this.nmf_fdg_home, ...
                    sprintf('baseline_%s_anticlust', opts.subgroup), ...
                    sprintf('baseline_%s_rep*', opts.subgroup)))';
                rec_errors = cell(1, length(pathDirs));
                loadData_ = @mladniArisCodes.loadData;
                calcRecError_ = @this.calcRecError;

                % parfor (rep = 1:length(pathDirs), opts.Ncores)
                for rep = 1:1
                    try



                        %% build the averaged image that NMF factors must estimate

                        Xmat_ = fullfile(pathDirs{rep}, sprintf('NumBases%i', param.numBase), 'X.mat');
                        if isfile(Xmat_) && opts.cache_files
                            load(Xmat_);  % while stored in, e.g., baseline_cn/NumBases24/X.mat, it contains averages
                            % pertinent to models of all cardinality
                        else

                            % ensure nifti_files_mounted.csv
                            inFiles_ = fullfile(pathDirs{rep}, 'nifti_files_mounted.csv');
                            if ~isfile(inFiles_)
                                g0 = strrep(inFiles_, '_mounted.csv', '.csv');
                                t = readtable(g0, Format='%s', Delimiter=' ', ReadVariableNames=false);
                                t.Var1 = strrep(t.Var1, '/scratch/jjlee', '/home/usr/jjlee');
                                writetable(t, inFiles_, WriteVariableNames=false);
                            end

                            % load imaging using Aris' conventions
                            data = loadData_(inFiles_, param, []);
                            meanX = mean(data.X, 2);
                            X = data.X; % data.X(meanX>0,:);
                            clear data;
                            save(Xmat_, 'X', 'meanX');
                        end
                        


                        %% collect rec_errors

                        resultsDir = fullfile(pathDirs{rep}, 'results');
                        ensuredir(resultsDir);                        
                        assert(size(X,1) == sum(selectVoxels), stackstr())
                        rec_errors{rep} = calcRecError_(X, pathDirs{rep}, resultsDir, selectVoxels);
                    catch ME
                        handwarning(ME)
                    end
                end
                save(fullfile(resultsDir_, [stackstr() '.mat']), 'rec_errors');
            else
                ld = load(fullfile(resultsDir_, [stackstr() '.mat']));
                rec_errors = ld.rec_errors;
            end





            
            % if isfile(this.Xmat) && this.cache_files
            %     load(this.Xmat);
            % else
            % 
            %     inFiles_ = fullfile(this.outputDir, 'nifti_files_mounted.csv');
            %     if ~isfile(inFiles_)
            %         g0 = strrep(inFiles_, '_mounted.csv', '.csv');
            %         t = readtable(g0, Format='%s', Delimiter=' ', ReadVariableNames=false);
            %         t.Var1 = strrep(t.Var1, '/scratch/jjlee', '/home/usr/jjlee');
            %         writetable(t, inFiles_, WriteVariableNames=false);
            %     end
            % 
            %     data = mladni.ArisCodes.loadData(inFiles_, param, []);
            %     meanX = mean(data.X,2);
            %     X = data.X; % data.X(meanX>0,:); 
            %     clear data;
            %     save(this.Xmat, 'X', 'meanX');
            % end
            % resultsDir = fullfile(this.outputDir, 'results');
            % ensuredir(resultsDir);
            % 
            % % BUG FIX:  selectVoxels must be logical, and
            % %           numel(selectVoxels) == numel(mask_for_imaging), and
            % %           size(selectVoxels) == [N_voxels 1], and
            % %           size(X,1) == sum(selectVoxels)
            % assert(isfile(param.mask))
            % mask_ic = mlfourd.ImagingContext2(param.mask);
            % selectVoxels = logical(mask_ic);
            % selectVoxels = selectVoxels(:);
            % assert(size(X,1) == sum(selectVoxels), stackstr())
            % this.calcRecError(X, this.outputDir, resultsDir, selectVoxels);
        end

        function [f,h] = build_rm_raincloud(this, opts)
            arguments
                this mlvg.Idif2024 
                opts.metric {mustBeTextScalar} = "metric"
                opts.fname {mustBeTextScalar} = ""
                opts.axis_label {mustBeTextScalar} = "axis_label"
                opts.idx_metric double = 1
                opts.rescale double = 1 
            end
            if isemptytext(opts.fname)
                opts.fname = lower(opts.metric);
            end

            %% colormapping

            [cb] = cbrewer2('Spectral', 12, 'pchip');
            cl(1, :) = cb(11, :);  % blue is IDIF
            cl(2, :) = cb(2, :);  % red is Twilite

            pwd0 = pushd(this.outputDir);

            %% aufbau data

            f = str2func(opts.fname);
            [nii_idif,mg_idif] = f(this, input_func="idif", stats="");
            [nii_twil,mg_twil] = f(this, input_func="twil", stats="");
            if ~isempty(getenv("DEBUG"))
                imagesc(nii_idif);
                imagesc(nii_twil);
            end

            % twilite data is more limited than idif data; use this.match_globbed()
            ifc_idif = nii_idif.imagingFormat;
            ifc_twil = nii_twil.imagingFormat;
            [ifc_idif, ifc_twil] = this.match_globbed(ifc_idif, mg_idif, ifc_twil, mg_twil);

            % select index of metric when necessary
            if ndims(ifc_idif) > 2 %#ok<ISMAT>
                ifc_idif.img = squeeze(ifc_idif.img(:,opts.idx_metric,:));
            end
            if ndims(ifc_twil) > 2 %#ok<ISMAT>
                ifc_twil.img = squeeze(ifc_twil.img(:,opts.idx_metric,:));
            end

            assert(all(size(ifc_idif.img) == size(ifc_twil.img)))
            N_sub = size(ifc_idif.img, ndims(ifc_idif.img));
            data = cell(N_sub, 2);  % N_sub x N. input func. methods
            for idx_sub = 1:N_sub
                data{idx_sub, 1} = opts.rescale * double(ifc_idif.img(:, idx_sub));
                data{idx_sub, 2} = opts.rescale * double(ifc_twil.img(:, idx_sub));
            end

            %% make & save figure

            f_metric  = figure(Position=this.fig_position);
            h = rm_raincloud(data, cl);
            %set(gca, 'YLim', [-0.3 1.6]);
            xlabel(opts.axis_label);
            if contains(opts.metric, ["OEF", "vcapillary", "fwatermetab", "CMRO2", "logZ_oo"])
                ylabel("Scans from 6 Participants");
            else
                ylabel("Participants");
            end
            fontsize(scale=2.5)
            saveFigure2(f_metric, fullfile(pwd, lower(opts.metric) + "_raincloud"));

            popd(pwd0);
        end
        
        function build_surfaces(this)
            pwd0 = pushd(fullfile(this.outputDir, sprintf('NumBases%i', this.selectedNumBases), 'OPNMF', 'niiImg'));
            mlshimony.Workbench.vol2surf('Basis_*.nii')
            popd(pwd0);
        end
        
        function T = build_table_variances(this, opts)
            arguments
                this mladni.NMF
                opts.subgroups {mustBeText}
                opts.numbases double = this.N_PATTERNS
            end
            subgroups = convertCharsToStrings(opts.subgroups);
            niidir = fullfile(this.nmf_fdg_home, "baseline_cn", "NumBases"+opts.numbases, "OPNMF", "niiImg");

            % NMF patterns
            for ib = 1:opts.numbases
                basis_ics(ib) = mlfourd.ImagingContext2(fullfile(niidir, "Basis_"+ib+".nii")); 
            end

            % table of variances for patterns (rows) and subgroups (cols)
            T = table();
            for isg = 1:length(subgroups)
                % T columns
                var_fns = glob(fullfile(this.nmf_fdg_home, subgroups(isg), "all_trc-FDG*pet_on_T1w_Warped_dlicv_variance.nii.gz"));
                assert(~isempty(var_fns))
                var_ic = mlfourd.ImagingContext2(var_fns{1});

                var_vec = nan(opts.numbases, 1);
                for ib = 1:opts.numbases
                    % T rows
                    var_vec(ib) = var_ic.volumeWeightedAveraged(basis_ics(ib)./dipsum(basis_ics(ib))); % weighted average
                end

                T = addvars(T, ascol(var_vec), NewVariableNames=subgroups(isg));
            end

            writetable(T, fullfile(this.nmf_fdg_home, stackstr()+".csv"));
            save(fullfile(this.nmf_fdg_home, stackstr()+".mat"), "T");
            heatmap(table2array(T));
            ylabel("NMF Patterns");
            xlabel("Diagnostic Groups");
            saveFigure2(gcf, stackstr());
        end
        
        function call(this)
            this.build_repeated_calc_rec_error()
            cd(this.outputDir);            
            this.evaluateRepeatedReproducibility2('cn', data_home=this.data_home)
        end
        
        function T = call2(this, tags)
            %% Writes imaging data to X2.mat, 
            %  then calls this.calculateSelectedComponentWeightedAverageNIFTI(), 
            %  which writes component_weighted_average_study-study_design.csv.

            arguments
                this mladni.NMF
                tags {mustBeTextScalar} = this.study_design
            end

            param.isList = 1 ;
            param.downSample = this.downSample;
            param.smooth = 0;
            param.mask = fullfile(this.data_home, 'VolBin/mask.nii.gz');
            param.permute = 0;
            param.numBase = this.selectedNumBases;
            param.componentDir = this.componentDir;

            assert(isfile(param.mask))

            if isfile(this.X2mat) && this.cache_files
                load(this.X2mat);
            else
                data = mladni.ArisCodes.loadData(this.inFiles2, param, []);
                meanX = mean(data.X,2);
                X = data.X;  % data.X(meanX>0,:); 
                clear data;
                save(this.X2mat, 'X', 'meanX' );
            end

            weightedAverFilename = fullfile(this.componentDir, sprintf('component_weighted_average_%s.csv', tags));
            this.calculateSelectedComponentWeightedAverageNIFTI( ...
                this.inFiles2, this.targetDatasetDir, this.selectedNumBases, weightedAverFilename);
            T = readtable(weightedAverFilename, Delimiter=',', ReadVariableNames=false);
        end
        
        function run_extractBasesMT(this, numBases, outputDir)
            assert(isscalar(numBases));
            assert(isfolder(outputDir));
            
            cmd = sprintf('%s/run_extractBasesMT.sh %s OPNMF %s 1 %i outputDir %s saveInterm 1 negPos 0 initMeth 4', ...
                this.volbin, this.mcrroot, this.inFiles, numBases, outputDir); % downSample 2
            system(cmd);
        end
        
        function run_nmf_dataset(this)
            %% quasi-constant

            this.mcrroot = '/data/nil-bluearc/raichle/jjlee/Local/MCR_zfs/R2018b/v95';
            if ~isempty(ipr.numBases)
                this.numBases = ipr.numBases;
            else
                for n = 2:2:mladni.NMF.MAX_NUM_BASES
                    if ~isfile(fullfile(this.outputDir, sprintf('NumBases%i', n), 'OPNMF', 'niiImg', ...
                            sprintf('Basis_%i.nii', n)))
                        try
                            cmd = sprintf('%s/run_nmf_dataset.sh -d %s -r %i', ...
                                this.volbin, this.nmfDataset, n);
                            fprintf(strcat('mladni.NMF: ', cmd))
                            mysystem(cmd);
                        catch ME
                            handerror(ME)
                        end                        
                    end
                    this.numBases = [this.numBases n];
                end
            end
        end
        
        function submit_nmf_dataset(this)
            cmd = sprintf('%s/submit_nmf_dataset.sh -d %s -m %i', this.volbin, this.nmfDataset, this.memInGB);
            [~,r] = mlbash(cmd);
            disp(r)
        end
    end

    methods (Static)
        function create_montage(varargin)
            %  Args:
            %      path (folder): e.g., .,
            %                          /path/to/nmf_dataset, 
            %                          /path/to/nmf_dataset/NumBases20, 
            %                          /path/to/nmf_dataset/NumBases20/OPNMF/niiImg, 
            
            import mladni.NMF;

            ip = inputParser;
            addOptional(ip, 'path', pwd, @isfolder);
            addParameter(ip, 'noclobber', false, @islogical)
            parse(ip, varargin{:});
            ipr = ip.Results;
            if strcmpi(mybasename(ipr.path), 'niiImg') % /path/to/niiImg
                pwd0 = pushd(ipr.path);
                atl = mlfourd.ImagingContext2( ...
                    fullfile(getenv('FSLDIR'), 'data', 'standard', 'MNI152_T1_1mm.nii.gz'));
                for g = glob('Basis_*.nii')'
                    try
                        ic = mlfourd.ImagingContext2(g{1});
                        product = strrep(strcat(ic.fqfp, '.png'), 'niiImg', 'Figures');
                        if ~isfile(product) || ...
                                (isfile(product) && ~ipr.noclobber)
                            atl.save_qc(ic);
                            png = strrep(g{1}, '.nii', '.png');
                            movefile(png, '../Figures', 'f');
                        end
                    catch ME
                        fprintf("%s:  %s\n", stackstr(), ME.message)
                    end
                end
                ensuredir('../Figures');
                popd(pwd0);
                return
            end
            if contains(basename(ipr.path), 'NumBases') % /path/to/NumBases*
                pth_to_niiImg = fullfile(ipr.path, 'OPNMF', 'niiImg', '');
                if isfolder(pth_to_niiImg)
                    NMF.create_montage(pth_to_niiImg, 'noclobber', ipr.noclobber);
                end
                return
            end
            if ~isempty(glob(fullfile(ipr.path, 'NumBases*'))) % nmf_dataset
                for g = glob(fullfile(ipr.path, 'NumBases*'))'
                    NMF.create_montage(g{1}, 'noclobber', ipr.noclobber);
                end
                return
            end
            warning('mladni:NotImplementedWarning', 'NMF.create_montage: nothing to be done');
        end
        
        function csv_out = filter_missing_images(varargin)
            %% reads csv, ignores missing rows, adjusts SINGULARITY_HOME as needed, writes csv_out
            
            ip = inputParser;
            addRequired(ip, 'csv', @isfile);
            addParameter(ip, 'csv_out', '', @istext)
            addParameter(ip, 'singularity_home', getenv('SINGULARITY_HOME'), @isfolder);
            parse(ip, varargin{:});
            ipr = ip.Results;  
            if isempty(ipr.csv_out)
                [pth,fp,x] = fileparts(ipr.csv);
                ipr.csv_out = fullfile(pth, strcat(fp, '_', x));
            end

            %fprintf('mladni.NMF.filter_missing_images():\n')
            
            lines = readlines(ipr.csv); % string vector
            lines1 = {};
            Nmissing = 0;
            for li = 1:length(lines)
                line_ = lines{li};
                if isempty(line_)
                    break
                end
                if ~contains(line_, ipr.singularity_home)
                    ss = strsplit(line_, 'Singularity');
                    line_ = strcat(ipr.singularity_home, ss{2});
                end
                if isfile(line_)
                    lines1 = [lines1; line_];
                else
                    fprintf('\tmissing %s\n', line_);
                    Nmissing = Nmissing + 1;
                end
            end

            writetable(table(lines1), ipr.csv_out, 'WriteVariableNames', false);
            if Nmissing > 0
                fprintf('\tdiscovered %i missing files\n', Nmissing);
            end
            fprintf('\tread %s\n', ipr.csv);
            fprintf('\twrote %s\n', ipr.csv_out);
            csv_out = ipr.csv_out;
        end
        
        function split_nifti_files(csv)
            %% DEPRECATED; see mladni.Adni.anticlust_cn_repeat

            [pth,fp,x] = myfileparts(csv);
            if isempty(pth)
                pth = pwd;
            end
            assert(isfolder(pth));
            assert(strcmp(x, '.csv'));
            pth1 = strcat(pth, '_split1');
            ensuredir(pth1)
            pth2 = strcat(pth, '_split2');
            ensuredir(pth2)

            % read table to split
            tbl = readtable(csv, 'ReadVariableNames', false, 'Delimiter', ' ');

            % define splits
            s1 = rand(size(tbl.Var1));
            s1 = s1 > 0.5;
            s2 = ~s1;

            % make splits
            Var1 = tbl.Var1(s1);
            csv1 = fullfile(pth1, strcat(fp, x));
            writetable(table(Var1), csv1, 'WriteVariableNames', false);

            Var2 = tbl.Var1(s2);
            csv2 = fullfile(pth2, strcat(fp, x));
            writetable(table(Var2), csv2, 'WriteVariableNames', false);
        end
        
        %% CHPC functions

        function propcluster()
            %% PROPCLUSTER.
            %  Currently, the script asks for 24Gb RAM, 4 CPUs and 24h of compute time. Asking for more CPUs will make 
            %  the executable finish faster. The more resources one asks, the more time it will probably take for the 
            %  job to start. If an executable needs more resources than the ones requested, it will fail. I am not sure 
            %  about the exact requirements here as this depends both on the number of samples and the dimensionality of 
            %  data. I am almost certain that the nmf computations will take more than 24h (especially for the higher 
            %  number of components). However, I would be sKeptical about requesting more time as it seems that the 
            %  queueing system heavily penalizes longer jobs (small number of max running jobs). Note that the code 
            %  saves intermediate results and can restart from these intermediate points by simply running the same 
            %  command and pointing to the same output directory.

            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '128000'; 
            c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = '';
            c.AdditionalProperties.WallTime = '24:00:00'; % 24 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        
        function parcall()
            %% DEPRECATED legacy method
            groups = mladni.NMF.groups;
            parfor gi = 1:length(groups)
                try
                    this = mladni.NMF( ...
                        nmfDataset=sprintf('baseline_%s', groups{gi}));
                    cd(this.outputDir);
                    call(this);
                    this.evaluateRepeatedReproducibility(groups{gi})
                catch ME
                    handwarning(ME)
                end
            end
        end
        
        function parcall2(opts)
            %% DEPRECATED legacy method
            
            arguments
                opts.data_home {mustBeFolder} = getenv('ADNI_HOME')
            end

            targetDatasetDir = fullfile(opts.data_home, 'NMF_FDG', 'baseline_cn');
            assert(isfolder(targetDatasetDir))

            b = 2:2:mladni.NMF.MAX_NUM_BASES;
            parfor bi = 1:length(b)
                groups = mladni.NMF.groups;
                for gi = 1:length(groups)
                    try
                        this = mladni.NMF( ...
                            selectedNumBases=b(bi), ...
                            nmfDataset=sprintf('baseline_%s', groups{gi}));
                        cd(this.outputDir);
                        call2(this, targetDatasetDir, 'on_cn');
                    catch ME
                        handwarning(ME)
                    end
                end
            end

            %% gather component files on local machine
            % for b in baseline_*; do
            % pushd $b
            % for n in NumBases*; do
            % pushd $n
            % rsync -ra login3.chpc.wustl.edu:/scratch/jjlee/Singularity/ADNI/NMF_FDG/$b/$n/components .
            % popd; done
            % popd; done
        end        
        
        function par_create_montage(gid)
            arguments
                gid {mustBeScalarOrEmpty} = 1
            end

            if gid
                groups = mladni.NMF.groups;
                this = mladni.NMF(nmfDataset=sprintf('baseline_%s', groups{gid}));
                cd(this.outputDir);
                this.create_montage(this.outputDir, noclobber=true);
            else
                groups = mladni.NMF.groups;
                parfor gi = 1:length(groups)
                    this = mladni.NMF(nmfDataset=sprintf('baseline_%s', groups{gi}));
                    cd(this.outputDir);
                    this.create_montage(this.outputDir, noclobber=true);
                end
            end
        end
        
        function getDebugLog(j,c)
            try
                c.getDebugLog(j)
            catch
                c.getDebugLog(j.Tasks(end))
            end
        end

        %% Aris' inference, called from other method functions
        
        function calculateComponentWeightedAverageNIFTI(dataList,resultsDir,numBases,outPath)
            %% loops through NumBases*, then accesses all /OPNMF/niiImg/Basis_*.nii to project files in dataList to Basis_*.nii
            % dataList: .csv file containing the images for which the coefficients need
            %           to be calculated. Full path for every image file is given in
            %           every line
            % resultsDir: directory where the NMF results are (i.e., the level where
            %             the NumBases folder is placed)
            % numBases: determines the solution for which one wants to calculate
            %           subject coefficients
            % outPath: fq-filenamne for saving output
            
            import mlniftitools.*;
            resultsDir = convertStringsToChars(resultsDir);
            outPath = convertStringsToChars(outPath);
            
            for b=1:length(numBases)

                fprintf("%s: numBases->%g\n", stackstr(), numBases(b))
                dataPath=[resultsDir '/NumBases' num2str(numBases) '/OPNMF/niiImg/'];

                % loading estimated non-negative components
                listing = glob(fullfile(dataPath, 'Basis_*.nii'));
                listing = listing(~contains(listing, 'argmax') & ~contains(listing, 'all'));

                if(length(listing)~=numBases(b))
                    error(['I cannot find ' num2str(numBases(b)) ' basis images in the folder ' dataPath ])
                end
                
                for i=1:numBases(b)
                    nii = load_untouch_nii([dataPath listing(i).name]);
                    B(:,i) = double(nii.img(:)');
                end
                
                % normalize to sum to 1
                Blen = sum(B,1);
                if any(Blen==0)
                    Blen(Blen==0) = 1;
                end
                nB = bsxfun(@times,B,1./Blen);
                
                % since the size and number of files is such that we can not manage
                % in batch mode, we are going to calculate weighted average values
                % subject by subject
                
                % read list
                fid=fopen(dataList,'r');
                if (fid == -1)
                    error(['extractBases:calculateComponentWeightedAverage ','Can not open ' list ' file.']);
                end
                datafullpath = textscan(fid,'%s\n');
                fclose(fid);
                
                datafullpath = datafullpath{1,1} ;
                datafullpath = cellstr(datafullpath) ;
                count = numel(datafullpath);
                
                fid = fopen(outPath,'w');
                frmtWrite='%s,';
                frmtWrite=[frmtWrite repmat('%f,',1,numBases(b)-1)]; %#ok<*AGROW>
                frmtWrite=[frmtWrite '%f\n'];
                
                wA = zeros(count,numBases) ;
                  
                %print statement added by tom
                disp(sprintf('calculateComponentWeightedAverageNIFTI: loading %d files', count)); %#ok<DSPS>
                for i=1:count                    
                    nii = load_untouch_nii(datafullpath{i});
                    wA(i,:) = double(nii.img(:)')*nB;
                    fprintf(fid,frmtWrite,datafullpath{i},wA(i,:)');
                end
                fclose(fid);
            end
        end
        
        function calculateSelectedComponentWeightedAverageNIFTI(dataLists,resultsDir,numBases,outFile)
            %% selects single folder NumBases\d, then accesses all /OPNMF/niiImg/Basis_*.nii to project files in dataList to Basis_*.nii
            % dataList: .csv file containing the images for which the coefficients need
            %           to be calculated. Full path for every image file is given in
            %           every line
            % resultsDir: directory where the NMF results are (i.e., the level where
            %             the NumBases folder is placed)
            % numBases: determines the solution for which one wants to calculate
            %           subject coefficients
            % outFile: CSV-file saving output
            
            import mlniftitools.*;
            resultsDir = convertStringsToChars(resultsDir);
            outFile = convertStringsToChars(outFile);
        
            % works with only single folder NumBases\d found from reproducibility analysis

            fprintf("%s: numBases->%g\n", stackstr(), numBases)
            dataPath=[resultsDir '/NumBases' num2str(numBases) '/OPNMF/niiImg/'];
            
            % loading estimated non-negative components
            listing = glob(fullfile(dataPath, 'Basis_*.nii'));
            listing = listing(~contains(listing, 'argmax') & ~contains(listing, 'all'));
            
            if(length(listing)~=numBases)
                error(['I cannot find ' num2str(numBases) ' basis images in the folder ' dataPath ])
            end
            
            for i=1:numBases
                nii = load_untouch_nii(listing{i});
                B(:,i) = double(nii.img(:)');
            end
            
            % normalize to sum to 1
            Blen = sum(B,1);
            if any(Blen==0)
                Blen(Blen==0) = 1;
            end
            nB = bsxfun(@times,B,1./Blen);
            
            % since the size and number of files is such that we can not manage
            % in batch mode, we are going to calculate weighted average values
            % subject by subject
            dataLists = ensureCell(dataLists);
            for dl = 1:length(dataLists)
    
                % read list
                fid=fopen(dataLists{dl},'r');
                if (fid == -1)
                    error(['extractBases:calculateComponentWeightedAverage ','Can not open ' list ' file.']);
                end
                datafullpath = textscan(fid,'%s\n');
                fclose(fid);
                
                datafullpath = datafullpath{1,1} ;
                datafullpath = cellstr(datafullpath) ;
                count = numel(datafullpath);
                
                fid = fopen(outFile,'w');
                frmtWrite='%s,';
                frmtWrite=[frmtWrite repmat('%f,',1,numBases-1)]; %#ok<*AGROW>
                frmtWrite=[frmtWrite '%f\n'];
                
                wA = zeros(count,numBases) ;

                for i=1:count                    
                    %print statement suggested by tom
                    if ~mod(i,10)
                        disp(sprintf('%s: (%d/%d)', stackstr(), i, count)); %#ok<DSPS>
                    end
                    nii = load_untouch_nii(datafullpath{i});
                    wA(i,:) = double(nii.img(:)')*nB;
                    fprintf(fid,frmtWrite,datafullpath{i},wA(i,:)');
                end
                fclose(fid);
            end
        end
        
        function [RecError,X] = calcBLSARecError(inFiles, resultDir)

            resultDir = convertStringsToChars(resultDir);

            % populating path
            % addpath(genpath('../../'));

            mfMeth = {'OPNMF','PCA','Jade'} ;

            % load data -> we need to estimate the age
            param.isList = 1 ;
            param.downSample = 1 ;
            data = mladni.ArisCodes.loadData(inFiles,param,[]) ;
            meanX = mean(data.X,2);
            X = data.X(meanX>0,:); clear data

            % hard coded results path
            listing = dir(resultDir);
            listing=listing(3:end) ;
            hh =cellfun(@(x) (strfind(x,'NumBases')==1),{listing(:).name},'UniformOutput',false) ;
            listing=listing(cellfun(@(x) (~isempty(x)&&(x==1)),hh));
            numDifBases=numel(listing) ;

            % sort them in ascending order
            basisNum = zeros(1,numDifBases) ;
            for i=1:numDifBases
                basisNum(i) = str2double(listing(i).name(9:end));    
            end
            [~,idx]=sort(basisNum) ;
            sortedBasisNum=basisNum(idx) ;

            % allocate memory 
            RecError=zeros(numDifBases,numel(mfMeth));

            for j=1:numel(mfMeth)
                for i=1:numDifBases
                    dataPath=[resultDir '/' listing(idx(i)).name '/' mfMeth{j} '/ResultsExtractBases.mat'];
                    % loading B (bases) and C (loading coefficients) - we do not
                    % need B
                    load(dataPath); %#ok<LOAD>
                    if(strcmp(mfMeth{j},'PCA'))
                        Est = B(meanX>0,:)*C + meanX(meanX>0)*ones(1,size(C,2)); clear B C
                    else
                        Est = B(meanX>0,:)*C ; clear B C
                    end
                    RecError(i,j) = norm(X-Est,'fro') ;
                    clear Est
                end
            end            
        end
        
        function RecError = calcRecError(X, resultDir, resultsDir, selectVoxels, opts)
            
            % Function that calculates the reconstruction error given data matrix X and the non-negative matrix 
            % factorization results saved in the resultDir directory. The function additionally saves figures that plot
            % the reconstruction error, the gradient of the reconstruction error, and the percentage of improvement as a 
            % function of the number of components
            
            % load results and calculate reconstruction error
            % We assume that the results directory is organized as follows: there is a folder for every solution (i.e.,
            % NumBases${K}, wehere K is the number of components for the solution), and inside the folder there is .mat 
            % file (ResultsExtractBases.mat) that contains the matrices W and H that were estimated by the non negative 
            % matrix factorization

            % selectVoxels is a logical vector of size N_voxels x 1, for N_voxels in the imaging field of view

            arguments
                X {mustBeNumeric}
                resultDir {mustBeFolder}
                resultsDir {mustBeTextScalar}
                selectVoxels logical
                opts.do_plot logical = false
            end

            resultDir = convertStringsToChars(resultDir); %% JJL
            resultsDir = convertStringsToChars(resultsDir); %% JJL    
            assert(islogical(selectVoxels))
            assert(size(X,1) == sum(selectVoxels, "all"), stackstr()) %% JJL's BUG CHECK
            
            listing = dir(resultDir);
            listing=listing(3:end);
            hh =cellfun(@(x) ( (strfind(x,'NumBases')==1)  ),{listing(:).name},'UniformOutput',false) ;
            listing=listing(cellfun(@(x) (~isempty(x)&&(x==1)),hh));
            numDifBases=numel(listing);
            
            % sort them in ascending order
            basisNum = zeros(1,numDifBases) ;
            for i=1:numDifBases
                basisNum(i) = str2double(listing(i).name(9:end));
            end
            [~,idx]=sort(basisNum) ;
            sortedBasisNum=basisNum(idx) ;
            
            RecError=zeros(numDifBases,1);
            for b=1:numDifBases
                disp(b/numDifBases)
                load( ...
                    fullfile(resultDir, ['NumBases', num2str(sortedBasisNum(b))], 'OPNMF', ...
                    'ResultsExtractBases.mat')) %#ok<LOAD>
                Est = B*C ;
                assert(size(Est,1) == size(selectVoxels,1), stackstr()) %% JJL's BUG CHECK
                Est = Est(selectVoxels,:); %% JJL
                RecError(b) = norm(X-Est,'fro') ;
                clear B C
            end
            
            if opts.do_plot
                plt = mladni.Plot();

                % make figures
                % 1) reconstruction error
                plt.plotxy(sortedBasisNum, RecError, ...
                    xlab="Number of patterns in model space", ...
                    ylab="Reconstruction error", ...
                    fileprefix=fullfile(resultsDir, "RecError"));

                % 2) gradient of reconstruction error
                plt.plotxy(sortedBasisNum(2:end), diff(RecError), ...
                    xlab="Number of patterns in model space", ...
                    ylab="Gradient of reconstruction error", ...
                    fileprefix=fullfile(resultsDir, "gradientRecError"));

                % 3) Fractional improvement over range of components used
                improvement = abs(RecError-RecError(1))./abs(RecError(1)-RecError(end));
                plt.plotxy(sortedBasisNum, improvement, ...
                    xlab="Number of patterns in model space", ...
                    ylab="Improvement with more patterns", ...
                    fileprefix=fullfile(resultsDir, "fractionalImprovementRecError"));

                %close all
            end
        end        
        
        function r = clustering_adjustedRand_fast(u,v)
            % clustering quality measures assumptions : 
            % 0 corresponds to background
            % we do not care about the background
            % cluster labels are assumed to be enumerated from 1 to max number of
            % clusters
            %
            % this function should not be used when comparing binary segmentations
            
            m=max(max(u),max(v));
            
            if(m == 1)
                error('ClusteringAdjustedRandFast:argChk','This method should not be used for comparing binary segmentations');
            end
            
            va=zeros(1,m);
            vb=zeros(1,m);
            mat=zeros(m);
            
            for i=1:m
                va(i) = sum(u==i) ;
                vb(i) = sum(v==i) ;
                hu = (u==i) ;
                for j=1:m
                    hv = (v==j) ;
                    mat(i,j) = sum(hu.*hv);
                end
            end
            
            ra=sum(sum(mat.*(mat-ones(m))))/2.0;
            
            rb=va*(va-ones(1,m))'/2.0;
            
            rc=vb*(vb-ones(1,m))'/2.0;
            
            rn=length(u)*(length(u)-1)/2.0;
            
            r=(ra-rb*rc/rn)/( 0.5*rb+0.5*rc-rb*rc/rn );
        end        
        
        function tf = components_are_available()
            tf = false;
        end        
        
        function evaluateRepeatedReproducibility(subgroup, N, K, opts)
            arguments
                subgroup {mustBeTextScalar} = 'cn'
                N double {mustBeInteger} = 50 % # repetitions
                K double {mustBeInteger} = 20 % # bases checked by VolBin/*
                opts.data_home {mustBeFolder} = getenv('ADNI_HOME')
                opts.do_plot logical = true
            end

            home = fullfile(opts.data_home, 'NMF_FDG');
            subgroups = { ...
                'cn', 'preclinical', ...
                'cdr_0p5_apos', ...
                'cdr_gt_0_aneg', 'cdr_gt_0p5_apos'};
            assert(any(strcmp(subgroup, subgroups)));
            
            outputDir_ = fullfile(home, sprintf('baseline_%s', subgroup));
            ensuredir(outputDir_);
            ARI = zeros(K, N);
            overlap = cell(K, N);
            for rep = 1:N
                try
                    pathDir1_ = fullfile(home, sprintf('baseline_%s_repA%i', subgroup, rep));
                    pathDir2_ = strrep(pathDir1_, 'repA', 'repB');
                    [~,~,A_,o_] = mladni.ArisCodes.evaluateReproducibility( ...
                        pathDir1_, pathDir2_, outputDir_, false);
                    ARI(:,rep) = A_;
                    for ik = 1:K
                        try
                            overlap{ik,rep} = o_{ik};
                        catch
                        end
                    end
                catch ME
                    handwarning(ME)
                end
            end

            mu_ARI = mean(ARI, 2)';
            sigma_ARI = std(ARI, 0, 2)';

            mu_overlap = cell(1, K);
            sigma_overlap = cell(1, K);
            for basis = 1:K
                try
                    accum = [overlap{basis, :}];
                    mu_overlap{basis} = mean(accum, "all");
                    sigma_overlap{basis} = std(accum, 0, "all");
                catch 
                end
            end
            mu_overlap = cell2mat(mu_overlap);
            sigma_overlap = cell2mat(sigma_overlap);
            save(fullfile(outputDir_, 'evaluateRepeatedReproducibility.mat'), ...
                'ARI', 'overlap', 'mu_ARI', 'sigma_ARI', 'mu_overlap', 'sigma_overlap');

            if opts.do_plot
                sortedBasisNum = 2:2:2*K;
                xconf = [sortedBasisNum sortedBasisNum(end:-1:1)];

                figure
                yconf = [mu_overlap+sigma_overlap mu_overlap(end:-1:1)-sigma_overlap];
                p = fill(xconf,yconf, 'red');
                p.FaceColor = [0.6 0.8 1];
                p.EdgeColor = 'none';
                hold on
                plot(sortedBasisNum, mu_overlap, 'bo-', 'LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',[0 0 0])
                hold off
                xlabel('Number of components','fontsize',12)
                ylabel({'Split-sample reproducibility';'mean +/- std of inner product'},'fontsize',12)
                title(strcat('baseline_', subgroup), Interpreter="none")
                set(gca,'fontsize',12)
                saveas(gcf, fullfile(outputDir_, 'innerProductReproducibility.fig'))
                saveas(gcf, fullfile(outputDir_, 'innerProductReproducibility.png'))
                
                figure
                yconf = [mu_ARI+sigma_ARI mu_ARI(end:-1:1)-sigma_ARI];
                p = fill(xconf,yconf, 'red');
                p.FaceColor = [0.6 0.8 1];
                p.EdgeColor = 'none';                
                hold on
                plot(sortedBasisNum, mu_ARI, 'bo-', 'LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',[0 0 0])
                hold off
                xlabel('Number of components','fontsize',12)
                ylabel({'Split-sample reproducibility' ;'mean +/- std of Adjusted Rand Index'},'fontsize',12)
                title(strcat('baseline_', subgroup), Interpreter="none")
                set(gca,'fontsize',12)
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility.fig'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility.png'))
            end
        end
        
        function evaluateRepeatedReproducibility2(subgroup, N, K, opts)
            arguments
                subgroup {mustBeTextScalar} = 'cn'
                N double {mustBeInteger} = 50 % # repetitions
                K double {mustBeInteger} = 20 % # bases checked by VolBin/*
                opts.fracx double = 0.5
                opts.Npx double = 3400
                opts.fracy double = 0.5
                opts.Npy double = 1440
                opts.data_home {mustBeFolder} = getenv('ADNI_HOME')
                opts.do_plot logical = true
            end

            import mladni.NMF.rm_raincloud

            nmf_fdg_home_ = fullfile(opts.data_home, 'NMF_FDG');
            subgroups = { ...
                'cn', 'preclinical', ...
                'cdr_0p5_apos', ...
                'cdr_gt_0_aneg', 'cdr_gt_0p5_apos'};
            assert(any(strcmp(subgroup, subgroups)));
            
            outputDir_ = fullfile(nmf_fdg_home_, sprintf('baseline_%s', subgroup), 'results');
            ensuredir(outputDir_);
            ARI = nan(K, N);
            overlap = cell(K, N);

            if ~isfile(fullfile(outputDir_, 'evaluateRepeatedReproducibility2.mat'))
                parfor rep = 1:N
                    try
                        pathDir1_ = fullfile(nmf_fdg_home_, sprintf('baseline_%s_repA%i', subgroup, rep));
                        pathDir2_ = strrep(pathDir1_, 'repA', 'repB');
                        [~,~,A_,o_] = mladni.ArisCodes.evaluateReproducibility( ...
                            pathDir1_, pathDir2_, outputDir_, false);
                        ARI(:,rep) = A_;
                        for ik = 1:K
                            try
                                overlap{ik,rep} = o_{ik};
                            catch %#ok<CTCH>
                            end
                        end
                    catch ME
                        handwarning(ME)
                    end
                end
                ARI = mladni.NMF.removeColsWithNans(ARI);
                overlap = mladni.NMF.removeColsWithNans(overlap);
                save(fullfile(outputDir_, 'evaluateRepeatedReproducibility2.mat'), ...
                    'ARI', 'overlap');
            else
                ld = load(fullfile(outputDir_, 'evaluateRepeatedReproducibility2.mat'));
                ARI = ld.ARI;
                overlap = ld.overlap;
            end

            if opts.do_plot
                sortedBasisNum = 2:2:2*K;
                sortedBasisNames = cellfun(@num2str, num2cell(2:2:2*K), UniformOutput=false);

                meanInner=cell2mat(cellfun(@(x) mean(x),overlap,'UniformOutput', false));
                meanInner(isnan(meanInner)) = 0;
                medianInner=cell2mat(cellfun(@(x) median(x),overlap,'UniformOutput', false));      
                medianInner(isnan(medianInner)) = 0;

                rm_raincloud(sortedBasisNum, sortedBasisNames, ARI');
                ylabel('Number of patterns in model space', 'fontsize', 30)
                xlabel({'Split-sample reproducibility'; '(Adjusted Rand Index)'}, 'fontsize', 30)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat2.fig'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat2.png'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat2.svg'))

                rm_raincloud(sortedBasisNum, sortedBasisNames, meanInner');
                ylabel('Number of patterns in model space', 'fontsize', 20)
                xlabel({'Split-sample reproducibility'; '(mean inner product)'}, 'fontsize', 20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat2.fig'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat2.png'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat2.svg'))

                rm_raincloud(sortedBasisNum, sortedBasisNames, medianInner');
                ylabel('Number of patterns in model space', 'fontsize', 20)
                xlabel({'Split-sample reproducibility'; '(median inner product)'}, 'fontsize', 20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat2.fig'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat2.png'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat2.svg'))

                ARI_snr = median(ARI', 1)./iqr(ARI', 1); %#ok<UDIM>
                plt = mladni.Plot();
                plt.plotxy(sortedBasisNum, ARI_snr, ...
                    xlab="Number of patterns in model space", ...
                    ylab=["Split-sample reproducibility"; "(adjusted Rand index, median/iqr)"], ...
                    fileprefix=fullfile(outputDir_, "MedianIqrARIReproducibility_repeat2"));
            end
        end
        
        function evaluateRepeatedReproducibility3(subgroup, N, K, opts)
            arguments
                subgroup {mustBeTextScalar} = 'cn'
                N double {mustBeInteger} = 50 % # repetitions
                K double {mustBeInteger} = 12 % # bases checked by VolBin/*
                opts.fracx double = 0.5
                opts.Npx double = 3400
                opts.fracy double = 0.5
                opts.Npy double = 1440
                opts.data_home {mustBeFolder} = fullfile(getenv('SINGULARITY_HOME'), 'OASIS3')
                opts.do_plot logical = true
            end

            import mladni.NMF.rm_raincloud

            nmf_fdg_home_ = fullfile(opts.data_home, 'NMF_FDG');
            
            outputDir_ = fullfile(nmf_fdg_home_, sprintf('baseline_%s', subgroup), 'results');
            ensuredir(outputDir_);
            ARI = nan(K, N);
            overlap = cell(K, N);

            if ~isfile(fullfile(outputDir_, 'evaluateRepeatedReproducibility3.mat'))
                oasis_folder = fullfile(nmf_fdg_home_, sprintf('baseline_%s', subgroup));
                anticlust_folders = glob(fullfile(nmf_fdg_home_, 'baseline_cn_repAdni*'));
                assert(length(anticlust_folders) >= N);
                parfor rep = 1:N
                    try
                        [~,~,A_,o_] = mladni.ArisCodes.evaluateReproducibility( ...
                            anticlust_folders{rep}, oasis_folder, outputDir_, false, K=K);
                        ARI(:,rep) = A_;
                        for ik = 1:K
                            try
                                overlap{ik,rep} = o_{ik};
                            catch %#ok<CTCH>
                            end
                        end
                    catch ME
                        handwarning(ME)
                    end
                end
                ARI = mladni.NMF.removeColsWithNans(ARI);
                overlap = mladni.NMF.removeColsWithNans(overlap);
                save(fullfile(outputDir_, 'evaluateRepeatedReproducibility3.mat'), ...
                    'ARI', 'overlap');
            else
                ld = load(fullfile(outputDir_, 'evaluateRepeatedReproducibility3.mat'));
                ARI = ld.ARI;
                overlap = ld.overlap;
            end

            if opts.do_plot
                sortedBasisNum = 2:2:2*K;
                sortedBasisNames = cellfun(@num2str, num2cell(2:2:2*K), UniformOutput=false);

                meanInner=cell2mat(cellfun(@(x) mean(x),overlap,'UniformOutput', false));
                meanInner(isnan(meanInner)) = 0;
                medianInner=cell2mat(cellfun(@(x) median(x),overlap,'UniformOutput', false));    
                medianInner(isnan(medianInner)) = 0;       

                rm_raincloud(sortedBasisNum, sortedBasisNames, ARI');
                ylabel('Number of patterns in model space', 'fontsize', 30)
                xlabel({'Split-sample reproducibility'; '(Adjusted Rand Index)'}, 'fontsize', 30)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat3.fig'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat3.png'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat3.svg'))

                rm_raincloud(sortedBasisNum, sortedBasisNames, meanInner');
                ylabel('Number of patterns in model space', 'fontsize', 20)
                xlabel({'Split-sample reproducibility'; '(mean inner product)'}, 'fontsize', 20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat3.fig'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat3.png'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat3.svg'))

                rm_raincloud(sortedBasisNum, sortedBasisNames, medianInner');
                ylabel('Number of patterns in model space', 'fontsize', 20)
                xlabel({'Split-sample reproducibility'; '(median inner product)'}, 'fontsize', 20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat3.fig'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat3.png'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat3.svg'))

                ARI_snr = median(ARI', 1)./iqr(ARI', 1); %#ok<UDIM>
                plt = mladni.Plot();
                plt.plotxy(sortedBasisNum, ARI_snr, ...
                    xlab="Number of patterns in model space", ...
                    ylab=["Split-sample reproducibility"; "(Adjusted Rand Index, median/iqr)"], ...
                    fileprefix=fullfile(outputDir_, "MedianIqrARIReproducibility_repeat3"));
            end
        end

        function rm_raincloud(pattern_nums, pattern_names, D, opts)
            %% pattern_nums is numeric.
            %  pattern_names is text.
            %  D is data ~ N_reps x K_patterns.

            arguments
                pattern_nums {mustBeNumeric}
                pattern_names {mustBeText}
                D {mustBeNumeric}
                opts.metric {mustBeTextScalar} = 'metric'
                opts.outputDir {mustBeFolder} = pwd
            end

            Nreps = size(D, 1);
            Kpatt = size(D, 2);
            fig_position = [80 80 1618 1000]; % coordinates for figures ~ [x0, y0, Dx, Dy]
            %cl = cbrewer2('Spectral', Kpatt, 'pchip');
            cl = [0.5, 0.5, 0.5];

            %% aufbau cell-array data

            data = cell(Kpatt, 1);
            for span = 1:Kpatt
                data{span, 1} = asrow(D(:, span));
            end

            %% make & save figure

            fig  = figure(Position=fig_position);
            h = rm_raincloud(data, cl);
            set(gca, 'YTickLabel', 2*Kpatt:-2:2);
            set(gca, 'XLim', [0, 1.1]);
            % xlabel(pattern_names); % defective
            % ylabel(opts.metric);
            fontsize(scale=2.5)
            % saveFigure2(fig, fullfile(opts.outputDir, lower(opts.metric) + "_raincloud"));
        end
        
        function data = loadDataFromList(list,param,subsetIdx)
            
            % Loading data form list
            % Input:
            %   list      : .txt file that contains function the full path for the data
            %               to be read
            %   param     : structure containing parameters
            %   subsetIdx : either EMPTY or FILENAME of the .txt file that contains
            %               that information regarding which images should be read
            %
            % This function returns a structure array data that has the following fields
            %   X         : matlab array that contains the actual data vectorized
            %               (dimension is D x N -> dimensionality x number of samples)
            %   y         : label information
            %   dimx      : dimension in x-axis
            %   dimy      : dimension in y-axis
            %   dimz      : dimension in z-axis
            %   nii       : nii structure of the input data that may be used in case
            %               one want to save the files in .nii format
            
            import mlniftitools.*;
            
            fid=fopen(list,'r');
            if (fid == -1)
                error(['extractBases:loadDataFromList ','Can not open ' list ' file.']);
            end
            datafullpath = textscan(fid,'%s %d\n');
            fclose(fid);
            
            data.y = datafullpath{1,2} ;
            datafullpath = datafullpath{1,1} ;
            datafullpath = cellstr(datafullpath) ;
            
            % keep filename info in the param structure, so that it can be saved
            data.filenames = datafullpath ;
            
            % find indexed entries, if trainIdx is not empty
            if(~isempty(subsetIdx))
                datafullpath = datafullpath(subsetIdx==1) ;
            end
            
            count = numel(datafullpath);
            info = load_untouch_header_only(datafullpath{1});
            
            data.dimx = info.dime.dim(2) ;
            data.dimy = info.dime.dim(3) ;
            data.dimz = info.dime.dim(4) ;
            
            % initialize data structure (vectorize data)
            if(param.downSample==1)
                % if no downsampling is performed
                disp('No downsampling is going to be performed to the data');
                % if data need to be masked, do not initialize yet the data structure
                if(isempty(param.mask))
                    data.X = zeros(data.dimx*data.dimy*data.dimz,count);
                end
                data.dsflag = 0 ;
            else
                % case of downsampling
                disp(['Data are going to be downsampled by a factor of ' num2str(param.downSample)]);
                new_size_x = ceil(data.dimx/param.downSample) ;
                new_size_y = ceil(data.dimy/param.downSample) ;
                new_size_z = ceil(data.dimz/param.downSample) ;
                % if data need to be masked, do not initialize yet the data structure
                if(isempty(param.mask))
                    data.X = zeros(new_size_x*new_size_y*new_size_z,count);
                end
                data.dsflag = 1;
            end
            
            % if data need to smoothed
            if(~isempty(param.smooth) && param.smooth ~= 0)
                disp('Data are going to be smoothed !');
                disp(['User specified smoothing kernel : ' num2str(param.smooth)]);
                [status,scratchDir]=tempdir; %system('echo ${SBIA_TMPDIR}');
                scratchDir=strtrim(scratchDir);
                if(status ~= 0)
                    error(['extractBases:loscratchDiradDataFromList ','Can not find which is the scratch directory.']);
                end
                if(isempty(scratchDir))
                    scratchDir='/tmp';
                end
                [status,pid]=system('echo $$');
                pid=str2double(pid);
                if(status ~= 0)
                    error(['extractBases:loadDataFromList ','Can not find process id.']);
                end
                % make temp directory to smooth images
                tmpDirName=[scratchDir '/smoothImg_NumBases_' num2str(param.numBase) '_PID_' num2str(pid) ] ;
                
                for tt=1:3
                    % if directory exists I will try to create a unique name 3 times
                    if(exist(tmpDirName,'dir'))
                        tmpDirName = [tmpDirName '_' num2str(floor(rand(1)*10000)) ] ;
                    end
                end
                if(~exist(tmpDirName,'dir'))
                    status=mkdir(tmpDirName);
                    if(status ~= 1)
                        error(['extractBases:loadDataFromList ','Can not create temporary directory ' tmpDirName ' to store smoothed images.']);
                    end
                else
                    error(['extractBases:loadDataFromList ','Someone else is writing in temporary directory ' tmpDirName ]);
                end
            else
                disp('No additional user specified smoothing is going to performed to the data');
            end
            
            % if data need to be masked
            if(~isempty(param.mask))
                disp('The data are going to be masked');
                disp(['Mask provided by user in ' param.mask]);
                mask_nii = load_untouch_nii(param.mask);
                mask = mask_nii.img ;
                
                % check that mask is given in the correct space
                if(mask_nii.hdr.dime.dim(2) ~= data.dimx || ...
                        mask_nii.hdr.dime.dim(3) ~= data.dimy || ...
                        mask_nii.hdr.dime.dim(4) ~= data.dimz )
                    error(['extractBases:loadDataFromList ','Mask size does not agree with image size !!']);
                end
                
                % initialize data structure
                if(param.downSample==1)
                    data.X = zeros(sum(mask(:)>0),count);
                else
                    xx1 = linspace(1,data.dimx,data.dimx) ;
                    yy1 = linspace(1,data.dimy,data.dimy) ;
                    zz1 = linspace(1,data.dimz,data.dimz) ;
                    xx2 = linspace(1,data.dimx,new_size_x) ;
                    yy2 = linspace(1,data.dimy,new_size_y) ;
                    zz2 = linspace(1,data.dimz,new_size_z) ;
                    [XX1,YY1,ZZ1] = meshgrid(yy1,xx1,zz1) ;
                    [XX2,YY2,ZZ2] = meshgrid(yy2,xx2,zz2) ;
                    mask = interp3(XX1,YY1,ZZ1,double(mask),XX2,YY2,ZZ2) ;
                    data.X = zeros(sum(mask(:)>0),count) ;
                end
            else
                disp('No user provided mask is specified.');
                mask = [] ;
            end
            data.mask = mask ;
            
            % load data
            % added print statement - Tom
            disp(sprintf("Loading %d images from %s ...", count, list)); %#ok<DSPS>
            for i=1:count                          
                try
                    if(~isempty(param.smooth) && param.smooth ~= 0)
                        smooth_command = ['3dmerge -1blur_fwhm ' num2str(param.smooth) ' -prefix ' tmpDirName '/smoothed_image_' num2str(i) '.nii.gz ' datafullpath{i} ] ;
                        [status,~] = system(smooth_command) ;
                        if(status ~= 0)
                            error(['extractBases:loadDataFromList ','Can not smooth image ' datafullpath{i} '.']);
                        end
                        % load smoothed image
                        nii = load_untouch_nii([tmpDirName '/smoothed_image_' num2str(i) '.nii.gz']);

                        % remove smoothed image
                        delete([tmpDirName '/smoothed_image_' num2str(i) '.nii.gz']);
                    else
                        nii = load_untouch_nii(datafullpath{i});
                    end
                    % if downsampling should be performed
                    if(param.downSample==1)
                        if(~isempty(mask))
                            data.X(:,i) = nii.img(mask>0) ;
                        else
                            data.X(:,i) = nii.img(:) ;
                        end
                    else
                        xx1 = linspace(1,size(nii.img,1),size(nii.img,1)) ;
                        yy1 = linspace(1,size(nii.img,2),size(nii.img,2)) ;
                        zz1 = linspace(1,size(nii.img,3),size(nii.img,3)) ;
                        xx2 = linspace(1,size(nii.img,1),new_size_x) ;
                        yy2 = linspace(1,size(nii.img,2),new_size_y) ;
                        zz2 = linspace(1,size(nii.img,3),new_size_z) ;
                        [XX1,YY1,ZZ1] = meshgrid(yy1,xx1,zz1) ;
                        [XX2,YY2,ZZ2] = meshgrid(yy2,xx2,zz2) ;
                        img_interp = interp3(XX1,YY1,ZZ1,double(nii.img),XX2,YY2,ZZ2) ;
                        if(~isempty(mask))
                            h_nii = double(img_interp(mask>0));
                            data.X(:,i) = h_nii(:) ; clear h_nii ;
                        else
                            data.X(:,i) = double(img_interp(:));
                        end
                    end
                catch ME
                    handwarning(ME)
                end
            end
            data.nii = nii ;
            
            % % if data were smoothed, remove directory
            if(~isempty(param.smooth) && param.smooth ~= 0)
                [status,~,~]=rmdir(tmpDirName);
                if(status ~= 1)
                    error(['mladni:loadDataFromList ','Can not remove folder ' tmpDirName ]);
                end
            end
        end
        
                
        function X = removeColsWithNans(X)
            if iscell(X)
                X = X(:, all(~isnan(cell2mat(X))));
                return
            end                        
            % remove reps (columns) that have any nan elements
            X = X(:, all(~isnan(X)));
        end
    end

    %% PRIVATE

    properties (Access = private)
        dataDemographics_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
