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
        fig_position = [80 80 2330 1440]; % coordinates for figures ~ [x0, y0, Dx, Dy]
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
            %  then repeatedly calls mladni.ArisCodes.calcRecError() for anticlustered bootstraps,
            %  then saves rec_errors for subsequent plotting.

            arguments
                this mladni.NMF
                opts.cache_files logical = true
                opts.subgroup {mustBeTextScalar} = 'cn'
                opts.N_cores {mustBeNumeric} = 46
            end

            param.isList = 1 ;
            param.downSample = this.downSample;
            param.smooth = 0;
            param.mask = char(fullfile(this.data_home, "VolBin", "mask.nii.gz"));
            param.permute = 0;
            param.numBase = this.selectedNumBases;
            param.componentDir = char(this.componentDir);

            % BUG FIX:  selectVoxels must be logical, and
            %           numel(selectVoxels) == numel(mask_for_imaging), and
            %           size(selectVoxels) == [N_voxels 1], and
            %           size(X,1) == sum(selectVoxels)
            assert(isfile(param.mask))
            mask_ic = mlfourd.ImagingContext2(param.mask);
            select_voxels = logical(mask_ic);
            select_voxels = select_voxels(:);

            results_dir_ = fullfile(this.nmf_fdg_home, sprintf("baseline_%s", opts.subgroup), "results");
            ensuredir(results_dir_);
            results_mat_ = fullfile(results_dir_, stackstr() + ".mat");

            %% build the estimates for each anticlustered bootstrap

            if ~isfile(results_mat_)
                path_dirs = mglob( ...
                    fullfile(this.nmf_fdg_home, ...
                    sprintf("baseline_%s_anticlust", opts.subgroup), ...
                    sprintf("baseline_%s_rep*", opts.subgroup)))';
                rec_errors = cell(1, length(path_dirs));
                loadData_ = @mladni.ArisCodes.loadData;
                calcRecError_ = @mladni.ArisCodes.calcRecError;

                progressbar('boostraps')
                % parfor (rep = 1:length(path_dirs), opts.N_cores)
                for rep = 1:length(path_dirs)
                    try

                        %% build the average of subject images which NMF patterns must estimate

                        X_mat_ = fullfile(path_dirs(rep), "NumBases" + param.numBase, "X.mat");
                        if isfile(X_mat_) && opts.cache_files
                            ld = load(X_mat_);  % while stored in, e.g., baseline_cn/NumBases24/X.mat, it contains averages
                                                % pertinent to models of all cardinality
                            X = ld.X;
                        else

                            % ensure nifti_files_mounted.csv
                            inFiles_ = fullfile(path_dirs(rep), "nifti_files_mounted.csv");
                            if ~isfile(inFiles_)
                                g0 = strrep(inFiles_, "_mounted.csv", ".csv");
                                t = readtable(g0, Format="%s", Delimiter=" ", ReadVariableNames=false);
                                t.Var1 = strrep(t.Var1, "/scratch/jjlee", "/home/usr/jjlee");
                                writetable(t, inFiles_, WriteVariableNames=false);
                            end

                            % load imaging using Aris' conventions
                            ldd = loadData_(inFiles_, param, []);
                            meanX = mean(ldd.X, 2);
                            X = ldd.X; % data.X(meanX>0,:);
                            % clear ldd;
                            save(X_mat_, "X", "meanX");
                        end

                        %% collect rec_errors

                        rep_results_dir_ = fullfile(path_dirs(rep), "results");
                        ensuredir(rep_results_dir_);                        
                        assert(size(X, 1) == sum(select_voxels), stackstr())
                        rec_errors{rep} = calcRecError_(X, path_dirs(rep), rep_results_dir_, select_voxels);

                        frep = (rep - 1) / length(path_dirs);
                        progressbar(frep)
                    catch ME
                        handwarning(ME)
                    end
                end
                save(results_mat_, "rec_errors");
            else
                ld = load(results_mat_);
                rec_errors = ld.rec_errors;
            end
        end

        function build_repeated_reproducibility(this, opts)
            %% assesses ADNI vs ADNI
            %
            %  Returns: struct S ~ 
            %     {ari_median, ari_iqr, ari_mean, ari_std, ...
            %      overlap_median, overlap_iqr, overlap_mean, overlap_std, ...
            %      }

            arguments
                this mladni.NMF
                opts.subgroup {mustBeTextScalar} = "cn"
                opts.N_rep double {mustBeInteger} = 50 % # bootstrapping repetitions
                opts.N_K double {mustBeInteger} = 20 % # models checked by VolBin/*
                opts.do_plot logical = true
            end

            results_dir = fullfile(this.nmf_fdg_home, "baseline_" + opts.subgroup, "results");
            ensuredir(results_dir);

            %% ARI, overlap_numer
            
            if isfile(fullfile(results_dir, stackstr() + ".mat"))

                % load
                ld = load(fullfile(results_dir, stackstr() + ".mat"));
                ARI = ld.ARI;
                overlap_numer = ld.overlap_numer;  % ~ 20 x 50
            else

                % build
                ARI = nan(opts.N_K, opts.N_rep);
                overlap = cell(opts.N_K, opts.N_rep);
                progressbar('boostraps', 'cardinality of models')
                for rep = 1:opts.N_rep
                    try
                        path_dir1 = fullfile( ...
                            this.nmf_fdg_home, ...
                            sprintf("baseline_%s_anticlust", opts.subgroup), ...
                            sprintf("baseline_%s_repA%i", opts.subgroup, rep));
                        path_dir2 = strrep(path_dir1, "repA", "repB");
                        [~,~,A_,o_] = mladni.ArisCodes.evaluateReproducibility( ...
                            path_dir1, path_dir2, results_dir, do_plot=false, rep=rep, N_rep=opts.N_rep);
                        ARI(:,rep) = A_;
                        for ik = 1:opts.N_K
                            try
                                overlap{ik,rep} = o_{ik};
                            catch %#ok<CTCH>
                            end
                        end
                    catch ME
                        handwarning(ME)
                    end
                end
                ARI = this.removeColsWithNans(ARI);
                overlap = this.removeColsWithNans(overlap);  % ~ 20 x 50 
                overlap_numer = cell2mat(cellfun(@median, overlap, UniformOutput=false));  % ~ 20 x 50
                overlap_numer(isnan(overlap_numer)) = 0;
                save(fullfile(results_dir, stackstr() + ".mat"), "ARI", "overlap_numer");
            end

            %% improvement of reconstruction errors

            if isfile(fullfile(results_dir, stackstr() + "_improvement.mat"))

                % load 
                ld = load(fullfile(results_dir, stackstr() + "_improvement.mat"));
                improvement_numer = ld.improvement_numer;
            else

                % build
                rec_errors_mat = fullfile(results_dir, "NMF_build_repeated_calc_rec_error.mat");
                assert(isfile(rec_errors_mat))
                ld = load(rec_errors_mat);
                rec_errors = ld.rec_errors;
                improvement = cell(1, opts.N_rep);
                for rep = 1:2*opts.N_rep
                    improvement{rep} = ...
                        abs(rec_errors{rep} - rec_errors{rep}(1)) ./ abs(rec_errors{rep}(1) - rec_errors{rep}(end));
                end
                improvement = this.removeColsWithNans(improvement);
                improvement_numer = cell2mat(improvement);  % N_K x N_rep
                save(fullfile(results_dir, stackstr() + "_improvement_numer.mat"), "improvement_numer");
            end

            %% do plot

            if opts.do_plot

                plt = mladni.Plot();
                span_model = 2:2:2*opts.N_K;
                d_span_model = 364.119 / 38;
                span_model_1 = (span_model - span_model(1)) * d_span_model;
                c_overlap = cbrewer2('YlOrRd', 8);
                c_ari = cbrewer2('YlGn', 5);
                c_recerr = cbrewer2('PuBu', 5);

                figure(Position=this.fig_position);
                hold on;

                % median overlap
                plt.rm_raincloud(overlap_numer', line_width=4, ...
                    rain_colours=c_overlap(end-3,:), pdf_colours=c_overlap(end-2,:), line_colour=c_overlap(end,:));

                % ARI
                plt.rm_plot(span_model_1, ARI', line_width=5, ...
                    colours=c_ari(end-1,:), line_colour=c_ari(end,:)); 

                % cumulative improvement of reconstruction error
                plt.rm_plot(span_model_1, improvement_numer', line_width=4, ...
                    colours=c_recerr(end-1,:), line_colour=c_recerr(end,:));                 

                hold off;
                xlabel("Reproducibility and reconstruction fidelity");
                ylabel("Number of patterns in model space");
                fontsize(scale=plt.fontsize_scale)
                text(0.805, 16*d_span_model, "adjusted Rand index of hard clusters", ...
                    FontSize=24, FontAngle="italic", Color=0.3*c_ari(end,:))
                text(0.41, 15*d_span_model, "inner product of matched patterns", ...
                    FontSize=24, FontAngle="italic", Color=0.3*c_overlap(end,:))
                text(0.1, 36*d_span_model, "cumulative improvement of reconstruction fidelity", ...
                    FontSize=24, FontAngle="italic", Color=0.3*c_recerr(end,:))
            end
        end
        
        function build_repeated_reproducibility2(this, subgroup, N_rep, N_K, opts)
            %% assesses OASIS3 vs OASIS3

            arguments
                this mladni.NMF
                subgroup {mustBeTextScalar} = 'cn'
                N_rep double {mustBeInteger} = 50 % # repetitions
                N_K double {mustBeInteger} = 20 % # bases checked by VolBin/*
                opts.fracx double = 0.5
                opts.Npx double = 3400
                opts.fracy double = 0.5
                opts.Npy double = 1440
                opts.data_home {mustBeFolder} = this.data_home
                opts.do_plot logical = true
            end

            nmf_fdg_home_ = fullfile(opts.data_home, 'NMF_FDG');
            subgroups = { ...
                'cn', 'preclinical', ...
                'cdr_0p5_apos', ...
                'cdr_gt_0_aneg', 'cdr_gt_0p5_apos'};
            assert(any(strcmp(subgroup, subgroups)));
            
            outputDir_ = fullfile(nmf_fdg_home_, sprintf('baseline_%s', subgroup), 'results');
            ensuredir(outputDir_);
            ARI = nan(N_K, N_rep);
            overlap = cell(N_K, N_rep);

            if ~isfile(fullfile(outputDir_, 'evaluateRepeatedReproducibility2.mat'))
                parfor rep = 1:N_rep
                    try
                        pathDir1_ = fullfile(nmf_fdg_home_, sprintf('baseline_%s_repA%i', subgroup, rep));
                        pathDir2_ = strrep(pathDir1_, 'repA', 'repB');
                        [~,~,A_,o_] = mladni.ArisCodes.evaluateReproducibility( ...
                            pathDir1_, pathDir2_, outputDir_, false);
                        ARI(:,rep) = A_;
                        for ik = 1:N_K
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
                plt = mladni.Plot();
                sortedBasisNum = 2:2:2*N_K;
                sortedBasisNames = cellfun(@num2str, num2cell(2:2:2*N_K), UniformOutput=false);

                meanInner=cell2mat(cellfun(@(x) mean(x),overlap,'UniformOutput', false));
                meanInner(isnan(meanInner)) = 0;
                medianInner=cell2mat(cellfun(@(x) median(x),overlap,'UniformOutput', false));      
                medianInner(isnan(medianInner)) = 0;

                plt.rm_raincloud(sortedBasisNum, sortedBasisNames, ARI');
                ylabel('Number of patterns in model space', 'fontsize', 30)
                xlabel({'Split-sample reproducibility'; '(Adjusted Rand Index)'}, 'fontsize', 30)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat2.fig'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat2.png'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat2.svg'))

                plt.rm_raincloud(sortedBasisNum, sortedBasisNames, meanInner');
                ylabel('Number of patterns in model space', 'fontsize', 20)
                xlabel({'Split-sample reproducibility'; '(mean inner product)'}, 'fontsize', 20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat2.fig'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat2.png'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat2.svg'))

                plt.rm_raincloud(sortedBasisNum, sortedBasisNames, medianInner');
                ylabel('Number of patterns in model space', 'fontsize', 20)
                xlabel({'Split-sample reproducibility'; '(median inner product)'}, 'fontsize', 20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat2.fig'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat2.png'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat2.svg'))

                ARI_snr = median(ARI', 1)./iqr(ARI', 1); %#ok<UDIM>
                plt.plotxy(sortedBasisNum, ARI_snr, ...
                    xlab="Number of patterns in model space", ...
                    ylab=["Split-sample reproducibility"; "(adjusted Rand index, median/iqr)"], ...
                    fileprefix=fullfile(outputDir_, "MedianIqrARIReproducibility_repeat2"));
            end
        end
        
        function build_repeated_reproducibility3(this, subgroup, N_rep, N_K, opts)
            %% assesses OASIS3 vs ADNI
            arguments
                this mladni.NMF
                subgroup {mustBeTextScalar} = 'cn'
                N_rep double {mustBeInteger} = 50 % # repetitions
                N_K double {mustBeInteger} = 12 % # bases checked by VolBin/*
                opts.fracx double = 0.5
                opts.Npx double = 3400
                opts.fracy double = 0.5
                opts.Npy double = 1440
                opts.data_home {mustBeFolder} = this.data_home
                opts.do_plot logical = true
            end

            nmf_fdg_home_ = fullfile(opts.data_home, 'NMF_FDG');
            
            outputDir_ = fullfile(nmf_fdg_home_, sprintf('baseline_%s', subgroup), 'results');
            ensuredir(outputDir_);
            ARI = nan(N_K, N_rep);
            overlap = cell(N_K, N_rep);

            if ~isfile(fullfile(outputDir_, 'evaluateRepeatedReproducibility3.mat'))
                oasis_folder = fullfile(nmf_fdg_home_, sprintf('baseline_%s', subgroup));
                anticlust_folders = glob(fullfile(nmf_fdg_home_, 'baseline_cn_repAdni*'));
                assert(length(anticlust_folders) >= N_rep);
                parfor rep = 1:N_rep
                    try
                        [~,~,A_,o_] = mladni.ArisCodes.evaluateReproducibility( ...
                            anticlust_folders{rep}, oasis_folder, outputDir_, false, K=N_K);
                        ARI(:,rep) = A_;
                        for ik = 1:N_K
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
                plt = mladni.Plot();                
                sortedBasisNum = 2:2:2*N_K;
                sortedBasisNames = cellfun(@num2str, num2cell(2:2:2*N_K), UniformOutput=false);

                meanInner=cell2mat(cellfun(@(x) mean(x),overlap,'UniformOutput', false));
                meanInner(isnan(meanInner)) = 0;
                medianInner=cell2mat(cellfun(@(x) median(x),overlap,'UniformOutput', false));    
                medianInner(isnan(medianInner)) = 0;       

                plt.rm_raincloud(sortedBasisNum, sortedBasisNames, ARI');
                ylabel('Number of patterns in model space', 'fontsize', 30)
                xlabel({'Split-sample reproducibility'; '(Adjusted Rand Index)'}, 'fontsize', 30)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat3.fig'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat3.png'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat3.svg'))

                plt.rm_raincloud(sortedBasisNum, sortedBasisNames, meanInner');
                ylabel('Number of patterns in model space', 'fontsize', 20)
                xlabel({'Split-sample reproducibility'; '(mean inner product)'}, 'fontsize', 20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat3.fig'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat3.png'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat3.svg'))

                plt.rm_raincloud(sortedBasisNum, sortedBasisNames, medianInner');
                ylabel('Number of patterns in model space', 'fontsize', 20)
                xlabel({'Split-sample reproducibility'; '(median inner product)'}, 'fontsize', 20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat3.fig'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat3.png'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat3.svg'))

                ARI_snr = median(ARI', 1)./iqr(ARI', 1); %#ok<UDIM>
                plt.plotxy(sortedBasisNum, ARI_snr, ...
                    xlab="Number of patterns in model space", ...
                    ylab=["Split-sample reproducibility"; "(Adjusted Rand Index, median/iqr)"], ...
                    fileprefix=fullfile(outputDir_, "MedianIqrARIReproducibility_repeat3"));
            end
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
            this.build_repeated_reproducibility('cn')
            % this.build_repeated_calc_rec_error()
        end
        
        function T = call2(this, tags)
            %% Writes imaging data to X2.mat, 
            %  then calls mladni.ArisCodes.calculateSelectedComponentWeightedAverageNIFTI(), 
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
            mladni.ArisCodes.calculateSelectedComponentWeightedAverageNIFTI( ...
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
                    this.build_repeated_reproducibility(groups{gi})
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
             
        function X = removeColsWithNans(X)
            %% implicitly converts cell -> mat

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
