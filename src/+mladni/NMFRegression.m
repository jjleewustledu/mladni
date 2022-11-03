 classdef NMFRegression < handle
    %% GAM Regressions with trees
    %  
    %  Created 14-Sep-2022 20:29:21 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2039608 (R2022a) Update 5 for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Constant)
        ARIS_COVARIATES = ...
            {'AmyloidStatus' 'MergeAge' 'MergePtGender' 'CDGLOBAL' ...
            'Components'}
        JACKS_COVARIATES = ...
            {'AmyloidStatus' 'MergeFdg' 'MergeApoE4' ...
            'MergeAge' 'MergePtGender' 'MergePtEducat' ...
            'CDGLOBAL' 'MergeCdrsbBl' 'MergeMmseBl' ...
            'Components'}
        SELECTED_COVARIATES = ...
            {'AmyloidStatus' 'MergeApoE4' ...
            'MergeAge' 'MergePtGender' 'MergePtEducat' ...
            'MergeCdrsbBl' 'MergeMmseBl' ...
            'CDMEMORY' 'CDORIENT' 'CDJUDGE' 'CDCOMMUN' 'CDHOME' 'CDCARE' 'CDGLOBAL' ...
            'Phase' 'SITEID' ...
            'Components'}
        MOST_COVARIATES = ...
            {'AmyloidStatus' 'MergeApoE4' 'ICV' 'PVE1' ...
            'MergeAge' 'MergePtGender' 'MergePtEducat' 'MergePtEthCat' 'MergePtRacCat' 'MergePtMarry' ...
            'MergeMmseBl' 'MergeCdrsbBl' ...
            'CDMEMORY' 'CDORIENT' 'CDJUDGE' 'CDCOMMUN' 'CDHOME' 'CDCARE' 'CDGLOBAL' ...
            'MergeAdas11Bl' 'MergeAdas13Bl' 'MergeAdasQ4Bl' ...
            'MergeRavltImmediateBl' 'MergeRavltLearningBl' 'MergeRavltForgettingBl' 'MergeRavltPercForgettingBl' ...
            'MergeEcogPtMemBl' 'MergeEcogPtLangBl' 'MergeEcogPtVisspatBl' 'MergeEcogPtPlanBl' 'MergeEcogPtOrganBl' ...
            'MergeEcogPtDivattBl' 'MergeEcogPtTotalBl' ...
            'MergeEcogSPMemBl' 'MergeEcogSPLangBl' 'MergeEcogSPVisspatBl' 'MergeEcogSPPlanBl' 'MergeEcogSPOrganBl' ...
            'MergeEcogSPDivattBl' 'MergeEcogSPTotalBl' ...
            'Phase' 'SITEID' ...
            'Components'}
    end

    properties
        componentDir
        do_bayesopt
        crossVal
        holdout
        is_bivariate
        predictors
        response
        reuse_mdl
        selected_covariates
        shape_functions
    end

    properties (Dependent)
        argmax_fileprefix
        component_weighted_average_csv
        idx_cn
        idx_dementia
        idx_mci
        Ncomp
        niiImgDir
        numBasesDir
        standardDir
        table1pp % for paper (table 1)++
        table_covariates
        table_covariates_mat
        table_features % table_selected w/o components, MergeCdrsbBl at end
        table_MergeCdrsbBl
        table_predictor
        table_selected % table_covariates w/ selected_covariates, response at end
        table_train
        table_test
    end

    methods

        %% GET

        function g = get.argmax_fileprefix(~)
            g = 'Basis_argmax_brain_mask_dil1';
        end
        function g = get.component_weighted_average_csv(this)
            g = fullfile(this.componentDir, 'component_weighted_average.csv');
        end        
        function g = get.idx_cn(this)
            g = strcmpi(this.table_train.MergeDx, 'CN');
        end
        function g = get.idx_dementia(this)
            g = strcmpi(this.table_train.MergeDx, 'Dementia');
        end
        function g = get.idx_mci(this)
            g = strcmpi(this.table_train.MergeDx, 'MCI');
        end
        function g = get.Ncomp(this)
            re = regexp(mybasename(this.numBasesDir), 'NumBases(?<N>\d+)', 'names');
            assert(~isempty(re))
            g = str2double(re.N);
        end
        function g = get.niiImgDir(this)
            g = fullfile(this.numBasesDir, 'OPNMF', 'niiImg');
        end
        function g = get.numBasesDir(this)
            g = fileparts(strtrim(this.componentDir));
        end
        function g = get.standardDir(~)
            g = fullfile(getenv('FSLDIR'), 'data', 'standard');
        end
        function g = get.table1pp(this)
            if ~isempty(this.table1pp_)
                g = this.table1pp_;
                return
            end
            table1pp_mat_ = fullfile(this.componentDir, 'table1pp.mat');
            if isfile(table1pp_mat_)
                ld = load(table1pp_mat_);
                this.table1pp_ = ld.table1pp;
                g = this.table1pp_;
                return
            end

            % build de novo
            pwd0 = pushd(this.componentDir);
            ld = load('pcor_cdrsb_and_components.mat');
            pcor = ld.pcor_cdrsb_and_components;
            ld = load('ms1model.mat'); % gam smooth
            ms1 = ld.ms1model;
            ld = load('ms3model.mat'); % linear
            ms3 = ld.ms3model;
            patterns = (1:this.Ncomp)';
            t = table(patterns);
            t.beta = ms3.beta;
            t.beta_se = ms3.se;
            t.lin_tval = ms3.tval;
            t.lin_pval = ms3.pval;
            t.lin_fdrp = ms3.fdrpval;
            t.edf = ms1.edf;
            t.rfDF = ms1.rfDF;
            t.F = ms1.F;
            t.smooth_pval = ms1.pval;
            t.smooth_fdrp = ms1.fdrpval;            
            t.pcor = pcor.pcor;
            t.pcor_pval = pcor.pval;
            t.pcor_tval = pcor.tval;
            t.pcor_gp = pcor.gp;
            t.n = pcor.n;            
            this.table1pp_ = t;
            popd(pwd0)
        end
        function g = get.table_covariates(this)
            if ~isempty(this.table_covariates_)
                g = this.table_covariates_;
                return
            end
            
            % find .mat
            if isfile(this.table_covariates_mat)
                ld = load(this.table_covariates_mat);
                this.table_covariates_ = ld.t;
                g = this.table_covariates_;
                return
            end
            % find .csv
            if isfile(this.component_weighted_average_csv)                
                fdg = mladni.FDGDemographics();
                this.table_covariates_ = fdg.table_covariates('component_weighted_average.csv');
                g = this.table_covariates_;
                return
            end
            error('mladni:IOError', ...
                'NMFRegression.get.table_covariates could not find dataframes')
        end
        function g = get.table_covariates_mat(this)
            g = fullfile(this.componentDir, 'mladni_FDGDemographics_table_covariates.mat');
        end
        function g = get.table_features(this)
            t = this.table_selected;
            t = this.removecomponents(t);
            g = movevars(t, 'MergeCdrsbBl', 'After', size(t,2));
        end
        function g = get.table_MergeCdrsbBl(this)
            t = this.table_covariates;

            w = table(t.Components, t.MergeCdrsbBl);
            w = renamevars(w, {'Var1' 'Var2'}, {'Component' 'MergeCdrsbBl'});
            g = splitvars(w);
        end
        function g = get.table_predictor(this)
            g = this.table1pp;
        end
        function g = get.table_selected(this)
            %% last column is response variable

            if ~isempty(this.table_selected_)
                g = this.table_selected_;
                return
            end

            g = table_optin(this, this.selected_covariates);
            g = movevars(g, this.response, 'After', size(g,2));
            %g = splitvars(g, 'Components');
        end
        function g = get.table_train(this)
            g = this.table_selected(training(this.cvpartition()),:);
        end
        function g = get.table_test(this)
            g = this.table_selected(test(this.cvpartition()),:);
        end

        %% 

        function ic = build_basis_argmax(this)
            %%

            pwd0 = pushd(this.niiImgDir);

            ic = mlfourd.ImagingContext2('Basis_1.nii');
            ic = ic.zeros;
            ifc = ic.nifti;
            for idx = 1:this.Ncomp
                ic_ = mlfourd.ImagingContext2(sprintf('Basis_%i.nii', idx));
                ic_ = ic_.blurred(3);
                ic_ = ic_.thresh(0.0001); % avoid argmax near eps
                ifc.img(:,:,:,idx) = ic_.nifti.img;
            end
            ifc.fileprefix = 'Basis_all'; % composite nifti of components
            ifc.save

            ic = mlfourd.ImagingContext2(ifc);
            ic_sumt = ic.timeSummed; % reconstruction from components
            ic_sumt = ic_sumt.blurred(2);
            ic_sumt = ic_sumt.thresh(0.0001);
            ic_sumt.save 
            
            ic_mask = ic_sumt.binarized; % mask of reconstruction from components
            ic_mask.save

            tmp = reshape(ifc.img, [91*109*91 this.Ncomp]);
            [~,argmax] = max(tmp'); %#ok<UDIM> 
            argmax = reshape(argmax, [91 109 91]);
            ifc_argmax = copy(ifc);
            ifc_argmax.img = argmax;
            ifc_argmax.fileprefix = 'Basis_argmax';
            ifc_argmax.save();

            % tight mask
            brain_mask = mlfourd.ImagingContext2(fullfile(this.standardDir, 'MNI152_T1_2mm_brain_mask.nii.gz'));
            ic_argmax = mlfourd.ImagingContext2(ifc_argmax);
            ic_argmax = ic_argmax .* brain_mask.binarized();
            ic_argmax.fileprefix = strcat(ifc_argmax.fileprefix, '_brain_mask');
            ic_argmax.save

            % dilated mask
            brain_mask = mlfourd.ImagingContext2(fullfile(this.standardDir, 'MNI152_T1_2mm_brain_mask_dil1.nii.gz'));
            ic_argmax = mlfourd.ImagingContext2(ifc_argmax);
            ic_argmax = ic_argmax .* brain_mask.binarized();
            ic_argmax.fileprefix = strcat(ifc_argmax.fileprefix, '_brain_mask_dil1');
            ic_argmax.save

            % more dilated mask
            brain_mask = mlfourd.ImagingContext2(fullfile(this.standardDir, 'MNI152_T1_2mm_brain_mask_dil.nii.gz'));
            ic_argmax = mlfourd.ImagingContext2(ifc_argmax);
            ic_argmax = ic_argmax .* brain_mask.binarized();
            ic_argmax.fileprefix = strcat(ifc_argmax.fileprefix, '_brain_mask_dil');
            ic_argmax.save

            popd(pwd0)
        end
        function ic = build_F_nifti(this, varargin)
            ip = inputParser;
            addParameter(ip, 'Ntoshow', 8, @isscalar);
            parse(ip, varargin{:});
            ipr = ip.Results;

            ifc = mlfourd.ImagingFormatContext2( ...
                fullfile(this.numBasesDir, 'OPNMF', 'niiImg', strcat(this.argmax_fileprefix, '.nii')));   
            img_ = zeros(size(ifc));
            t = sortrows(this.table1pp, "F", "descend");
            for c = 1:ipr.Ntoshow
                img_(ifc.img == t.patterns(c)) = t.F(c);
            end
            ifc.img = img_;
            ifc.fqfilename = fullfile(this.componentDir, 'F.nii.gz');
            save(ifc);
            ic = mlfourd.ImagingContext2(ifc);
        end
        function t = print_table1(this)
            t = removevars(this.table1pp, ...
                {'lin_pval' 'lin_fdrp' 'rfDF' 'smooth_pval' 'smooth_fdrp' 'pcor_pval' 'pcor_tval' 'pcor_gp'});

            % remaining:  patterns, beta, beta_se, lin_tval, edf, F, pcor
            for r = 1:size(t,1)
                fprintf('%i & %0.2g & %0.2g & %0.3g & %0.3g & %0.3g & %0.2g \\\\ \n', ...
                    t.patterns(r), t.beta(r), t.beta_se(r), t.lin_tval(r), t.edf(r), t.F(r), t.pcor(r));
            end
        end
        function ic = build_tval_nifti(this, varargin)
            ip = inputParser;
            addParameter(ip, 'Ntoshow', this.Ncomp, @isscalar);
            parse(ip, varargin{:});
            ipr = ip.Results;

            ifc = mlfourd.ImagingFormatContext2( ...
                fullfile(this.numBasesDir, 'OPNMF', 'niiImg', strcat(this.argmax_fileprefix, '.nii')));            
            img_ = zeros(size(ifc));
            t = sortrows(this.table1pp, "lin_tval", "ascend");
            for c = 1:ipr.Ntoshow
                img_(ifc.img == t.patterns(c)) = t.lin_tval(c);
            end
            ifc.img = img_;
            if ipr.Ntoshow < this.Ncomp
                fn = sprintf('lin_tval_best%i.nii.gz', ipr.Ntoshow);
            else
                fn = 'lin_tval.nii.gz';
            end
            ifc.fqfilename = fullfile(this.componentDir, fn);
            save(ifc);
            ic = mlfourd.ImagingContext2(ifc);
        end
        function cmdl = compact(this)
            if ~isempty(this.compact_mdl_)
                cmdl = this.compact_mdl_;
                return
            end
            this.compact_mdl_ = compact(this.mdl);
            cmdl = this.compact_mdl_;
        end
        function cvp = cvpartition(this)
            %% cross-validation partitions

            if ~isempty(this.cvp_)
                cvp = this.cvp_;
                return
            end

            matfile = strcat(clientname(true, 2), '.mat');
            if isfile(matfile) && this.reuse_mdl
                ld = load(matfile);
                this.cvp_ = ld.cvp;
                cvp = this.cvp_;
                return
            end

            this.cvp_ = cvpartition(size(this.table_selected, 1), 'HoldOut', this.holdout);
            cvp = this.cvp_;
            save(matfile, 'cvp');
        end
        function mdlc = fitlm(this, indices, modelspec, vars)
            arguments
                this mladni.NMFRegression
                indices {mustBeInteger} = 5
                modelspec {mustBeTextScalar} = 'quadratic'
                vars cell {mustBeText} = this.ARIS_COVARIATES
            end
            
            mdlc = cell(size(indices));
            t = table_optin(this, vars);
            u = this.removecomponents(t);          
            for idx = indices
                Pattern = t.Components(:,idx); 
                v_ = addvars(u, Pattern); % added to far right of table
                mdl_ = fitlm(v_, modelspec, PredictorVar=vars, ResponseVar='Pattern');
                mdlc{idx} = mdl_; 
                fprintf('Pattern %i  ============================================================================================================================================================\n', idx)
                Coefficients = mdl_.Coefficients;
                Coefficients(Coefficients.pValue > 0.05/this.Ncomp,:) = [];
                disp(sortrows(Coefficients, 'pValue'))
            end
            nonempty = cellfun(@(x) ~isempty(x), mdlc);
            mdlc = mdlc(nonempty);
        end
        function mdlc = fitlm_ArisDx(this, indices)
            predictors_ = {'MergeAge' 'MergePtGender' 'ArisDx'};

            mdlc = cell(size(indices));
            t = this.table_covariates;
            Components = t.Components;
            t = mladni.NMFRegression.removenuisance(t);
            u = mladni.NMFRegression.removecomponents(t);
            
            ArisDx = cell(size(u.MergePtGender));
            select = u.AmyloidStatus == 0 & u.CDGLOBAL == 0;
            ArisDx(select) = {'cn'};
            select = u.AmyloidStatus == 1 & u.CDGLOBAL == 0;
            ArisDx(select) = {'preclinical'};
            select = u.AmyloidStatus == 1 & u.CDGLOBAL == 0.5;
            ArisDx(select) = {'cdr0p5'};
            select = u.AmyloidStatus == 1 & u.CDGLOBAL > 0.5;
            ArisDx(select) = {'cdr1pp'};
            u = addvars(u, ArisDx);
            totoss = cellfun(@isempty, ArisDx);
            u(totoss,:) = [];
            Components(totoss,:) = [];
            
            for idx = indices
                Pattern = Components(:,idx); 
                v_ = addvars(u, Pattern); % added to far right of table
                mdl_ = fitlm(v_, 'quadratic', PredictorVar=predictors_, ResponseVar='Pattern');
                mdlc{idx} = mdl_; 
                fprintf('Pattern %i  ==================================================================================================\n', idx)
                Coefficients = mdl_.Coefficients;
                Coefficients(Coefficients.pValue > 0.05/this.Ncomp,:) = [];
                disp(sortrows(Coefficients, 'pValue'))
            end
            nonempty = cellfun(@(x) ~isempty(x), mdlc);
            mdlc = mdlc(nonempty);
        end
        function [mdlc,explainer] = fitrgam(this, indices, opts)
            %  Returns:
            %      mdlc RegressionGAM

            arguments
                this mladni.NMFRegression
                indices {mustBeInteger} = 5
                opts.cdr double = 0.5
                opts.queryPoint table = []
                opts.savefigs logical = false
                opts.close_fig logical = false
            end

            mdlc = cell(size(indices));
            explainer = cell(size(indices));
            t = this.table_train;
            u = mladni.NMFRegression.removecomponents(t);            
            for idx = indices
                Pattern = t.Components(:,idx); 
                v_ = addvars(u, Pattern); % added to far right of table
                mdl_ = fitrgam(v_, 'Pattern', 'FitStandardDeviation', true, ...
                    'OptimizeHyperparameters', 'all-univariate', ...
                    'HyperparameterOptimizationOptions', ...
                    struct('AcquisitionFunctionName','expected-improvement-plus','ShowPlots',false,'Verbose',0));
                mdlc{idx} = mdl_; 
                disp(mdl_)
                if isempty(opts.queryPoint)
                    opts.queryPoint = this.queryPoint(v_, idx, opts.cdr);
                end
                explainer{idx} = shapley(mdl_, QueryPoint=opts.queryPoint, UseParallel=true);
                disp(explainer{idx}.ShapleyValues)
                figure('DefaultAxesFontSize',14)
                N = length(this.selected_covariates);
                for nvars = N:-1:1
                    try
                        plot(explainer{idx}, NumImportantPredictors=nvars)
                        break
                    catch
                    end
                end
                xlabel('Shapley Value', FontWeight='bold');
                ylabel('Predictor', FontWeight='bold');
                if opts.cdr < 1
                    ti = sprintf('Shapley Explanation for Pattern %i; query CDR %g', idx, opts.cdr);
                else
                    ti = sprintf('Shapley Explanation for Pattern %i; query CDR 1:2', idx);
                end
                title(ti)
                %view(mdl_, "Mode", "graph")
            end
            nonempty = cellfun(@(x) ~isempty(x), mdlc);
            mdlc = mdlc(nonempty);

            if opts.savefigs
                p = sprintf('fitrgam_shapley_cdr%g_pattern', opts.cdr);
                p = strrep(p, '.', 'p');
                saveFigures(pwd, 'closeFigure', opts.close_fig, 'prefix', p)

                save(sprintf('fitrgam_mdlc_cdr%g.mat', opts.cdr), 'mdlc');
                save(sprintf('fitrgam_explainer_cdr%g.mat', opts.cdr), 'explainer');
            end
        end
        function mdlc = fitrtree(this, indices, opts)
            %  Returns:
            %      mdlc RegressionTree

            arguments
                this mladni.NMFRegression
                indices {mustBeInteger} = 5
                opts.cdr double = 0.5
                opts.queryPoint table = []
                opts.savefigs logical = false
                opts.close_fig logical = false
            end

            mdlc = cell(size(indices));
            t = this.table_train;
            u = mladni.NMFRegression.removecomponents(t);            
            for idx = indices
                Pattern = t.Components(:,idx); 
                v_ = addvars(u, Pattern); % added to far right of table
                mdl_ = fitrtree(v_, 'Pattern');
                mdlc{idx} = mdl_; 
                disp(mdl_)
                if isempty(opts.queryPoint)
                    opts.queryPoint = this.queryPoint(v_, idx, opts.cdr);
                end
                explainer = shapley(mdl_, QueryPoint=opts.queryPoint, UseParallel=true);
                disp(explainer.ShapleyValues)
                figure('DefaultAxesFontSize',14)
                N = length(this.selected_covariates);
                for nvars = N:-1:1
                    try
                        plot(explainer, NumImportantPredictors=nvars)
                        break
                    catch
                    end
                end
                xlabel('Shapley Value', FontWeight='bold');
                ylabel('Predictor', FontWeight='bold');
                if opts.cdr < 1
                    ti = sprintf('Shapley Explanation for Pattern %i; query CDR %g', idx, opts.cdr);
                else
                    ti = sprintf('Shapley Explanation for Pattern %i; query CDR 1:2', idx);
                end
                title(ti)
                %view(mdl_, "Mode", "graph")
            end
            nonempty = cellfun(@(x) ~isempty(x), mdlc);
            mdlc = mdlc(nonempty);

            if opts.savefigs
                p = sprintf('fitrtree_shapley_cdr%g_pattern', opts.cdr);
                p = strrep(p, '.', 'p');
                saveFigures(pwd, 'closeFigure', opts.close_fig, 'prefix', p)
            end
        end
        function gam = generalized_additive_model(this, varargin)
            %% web(fullfile(docroot, 'stats/regressiongam.html#mw_6b1cf57e-6602-4aa4-a76a-ee199299b0b7'))
            %  Params:
            %      response (any):
            %      predictors (any):
            %      shape_functions (any):
            %  Returns:
            %      gam: RegressionGAM object

            ip = inputParser;
            addParameter(ip, 'response', []);
            addParameter(ip, 'predictors', []);
            addParameter(ip, 'shape_functions', []);
            parse(ip, varargin{:})
            ipr = ip.Results;


        end        
        function m = mdl(this)
            if this.do_bayesopt
                m = this.mdl_bayesopt();
                return
            end
            if strcmpi(this.crossVal, 'on')
                m = this.mdl_cv();
                return
            end
            if this.is_bivariate
                m = this.mdl_bivariate();
            else
                m = this.mdl_univariate();
            end
        end
        function histogramPartialDependence(this)
            %% web(fullfile(docroot, 'stats/regressiontree.partialdependence.html'))

            tic
            CMdl = compact(this);
            fprintf([clientname(true, 2) ':  ']);
            toc

            numPoints = 20;
            Age = this.table_train.MergeAge;            
            ptX = linspace(min(Age), max(Age),numPoints)';

            for p = CMdl.PredictorNames(3:end)
                Component = this.table_train.(p{1});
                ptY = linspace(min(Component), max(Component),numPoints)';
                [pd,x,y] = partialDependence(CMdl, ["MergeAge", p{1}], this.table_train, "QueryPoints", [ptX ptY]);

                figure;
                t = tiledlayout(5,5,"TileSpacing","compact");

                ax1 = nexttile(2,[4,4]);
                imagesc(x,y,pd)
                title(strcat("Partial Dependence Plot for ", p{1}), "Interpreter", "none")
                colorbar("eastoutside")
                ax1.YDir = "normal";
                
                ax2 = nexttile(22,[1,4]);
                dX = diff(ptX(1:2));
                edgeX = [ptX-dX/2;ptX(end)+dX];
                histogram(Age,edgeX);
                xlabel("Age")
                xlim(ax1.XLim);
                
                ax3 = nexttile(1,[4,1]);
                dY = diff(ptY(1:2));
                edgeY = [ptY-dY/2;ptY(end)+dY];
                histogram(Component,edgeY)
                xlabel(p{1}, "Interpreter", "none")
                xlim(ax1.YLim);
                ax3.XDir = "reverse";
                camroll(-90)
            end
 
            saveFigures(pwd, 'closeFigure', false, 'prefix', 'partial dependence histogram Components_')
        end
        function h = plot(this, varargin)
            %% web(fullfile(docroot, 'stats/regressiongam.predict.html')) -> "Plot Prediction Intervals"
            %
            %  Params:
            %      X (table|numeric):  predictor data associated with GAM model

            ip = inputParser;
            ip.KeepUnmatched = true;
            addOptional(ip, 'X', this.table_test, @(x) istable(x) || isnumeric(x));
            addParameter(ip, 'Alpha', 0.01, @isscalar);
            parse(ip, varargin{:});
            ipr = ip.Results;

            [yFit,~,yInt] = predict(this, ipr.X, 'Alpha', ipr.Alpha, varargin{:});
            
            figure
            yTrue = ipr.X.(this.response);
            [sortedYTrue,I] = sort(yTrue); 
            plot(sortedYTrue,'o')
            hold on
            plot(yFit(I))
            try
                plot(yInt(I,1),'k:')
                plot(yInt(I,2),'k:')
                legend('True responses','Predicted responses', ...
                    'Prediction interval limits','Location','best')
            catch
                legend('True responses','Predicted responses')
            end
            ylabel(this.response)
            xlabel('subject')
            hold off
        end
        function ple = plotLocalEffects(this, varargin)
            %  Params:
            %      arow (scalar):  default := 1.
            %      new_figure (logical):  default := true.
            %      num_terms (scalar):  default := length(CMdl.PredictorNames).
            %  Returns:
            %      ple: bar object, web(fullfile(docroot, 'matlab/ref/matlab.graphics.chart.primitive.bar-properties.html?browser=F1help'))

            CMdl = compact(this);

            ip = inputParser;
            addParameter(ip, "arow", 1, @(x) 0 < x && x < size(this.table_train,1) + 1)
            addParameter(ip, "new_figure", true, @islogical)
            addParameter(ip, "num_terms", length(CMdl.PredictorNames), @isscalar)
            parse(ip, varargin{:});
            ipr = ip.Results;

            if ipr.new_figure
                figure
            end
            ple = plotLocalEffects(CMdl,this.table_train(ipr.arow,:), ...
                'IncludeIntercept', true, ...
                'NumTerms', ipr.num_terms);
            ple.Parent.TickLabelInterpreter = 'none';
            saveFigures(pwd, 'closeFigure', false, 'prefix', clientname(true, 2))
        end
        function plotPartialDependence(this)
            %% web(fullfile(docroot, 'stats/regressiontree.plotpartialdependence.html?browser=F1help'))

            tic
            CMdl = compact(this);
            fprintf([clientname(true, 2) ':  ']);
            toc

            % Create a partial dependence plot for the terms MergeAge and {Components_*} using the train & test data.
            figure
            plotPartialDependence(CMdl,["MergeAge","MergePtGender"], this.table_train)
            for p = CMdl.PredictorNames(2:end)
                tic
                figure
                plotPartialDependence(CMdl,["MergeAge", p{1}], this.table_train)
                toc
            end
        end
        function [yFit,ySD,yInt] = predict(this, X, varargin)
            %% web(fullfile(docroot, 'stats/regressiongam.predict.html'))
            %
            %  Params:
            %      X (table|numeric):  predictor data associated with GAM model
            %  Returns:
            %      vector of predicted response for predictor data in X.

            ip = inputParser;
            ip.KeepUnmatched = true;
            addOptional(ip, "X", this.table_test, @(x) istable(x) || isnumeric(x));
            parse(ip, X, varargin{:});
            ipr = ip.Results;

            if isa(this.mdl, 'RegressionPartitionedGAM')
                yFit = kfoldPredict(this, varargin{:});
                ySD = [];
                yInt = [];
                return
            end

            cmdl = compact(this.mdl);
            try
                [yFit,ySD,yInt] = predict(cmdl, ipr.X, varargin{:});
            catch
                yFit = predict(cmdl, ipr.X, varargin{:});
                ySD = [];
                yInt = [];
            end
        end
        function predictorImportance(this, tree)
            %% web(fullfile(docroot, 'stats/compactregressiontree.predictorimportance.html?browser=F1help'))

            arguments
                this mladni.NMFRegression %#ok<INUSA> 
                tree RegressionTree
            end
            
            imp = predictorImportance(tree);
            figure;
            bar(imp);
            title(sprintf('Fitted Regression Tree for %s', tree.ResponseName));
            ylabel('Predictor Importance Estimates', 'FontWeight', 'bold');
            xlabel('Predictors', 'FontWeight', 'bold');
            h = gca;
            h.XTickLabel = tree.PredictorNames;
            h.TickLabelInterpreter = 'none';
            h.XTickLabelRotation = 90;
            len = length(tree.PredictorNames);
            set(gca, 'XTick', 1:len, 'XTickLabel', tree.PredictorNames)
        end
        function this = regressionLearner(this)
            %% web(fullfile(docroot, 'stats/regressionlearner-app.html'))

            regressionLearner();  
        end        
        function mdlc = stepwiselm(this, indices, modelspec, vars)
            arguments
                this mladni.NMFRegression
                indices {mustBeInteger} = 5
                modelspec {mustBeTextScalar} = 'quadratic'
                vars cell {mustBeText} = this.JACKS_COVARIATES
            end

            mdlc = cell(size(indices));
            t = table_optin(this, vars);
            u = mladni.NMFRegression.removecomponents(t);
            for idx = indices
                Pattern = t.Components(:,idx); 
                v_ = addvars(u, Pattern); % added to far right of table
                mdl_ = stepwiselm(v_, modelspec, ResponseVar='Pattern', Criterion='bic');
                mdlc{idx} = mdl_; 
                fprintf('Pattern %i  ============================================================================================================================================================\n', idx)
                Coefficients = mdl_.Coefficients;
                Coefficients(Coefficients.pValue > 0.05/this.Ncomp,:) = [];
                disp(sortrows(Coefficients, 'pValue'))
            end
            nonempty = cellfun(@(x) ~isempty(x), mdlc);
            mdlc = mdlc(nonempty);
        end
        function t = table_optin(this, vars)
            arguments
                this mladni.NMFRegression
                vars cell {mustBeText} = this.selected_covariates
            end

            t = table;
            t_varnames = {};
            t_ = this.table_covariates;
            t_ = this.removenuisance(t_);
            t__varnames = t_.Properties.VariableNames;
            for v = t__varnames
                if contains(v{1}, vars)
                    t = addvars(t, t_.(v{1}));
                    t_varnames = [t_varnames v{1}]; %#ok<AGROW> 
                end
            end
            t.Properties.VariableNames = t_varnames;
        end

        function this = NMFRegression(varargin)
            %% NMFREGRESSION 
            %  Args:
            %      componentDir (folder):  containing component_weighted_average.csv
            %                              generatead by mladni.NMF.call(), 
            %                              e.g., baseline4/NumBases32/components
            %      crossVal (text):  default := 'on'.
            %      do_bayesopt (logical):  default := true
            %      holdout (scalar):  < 1 for training, testing with cross-validation.
            %      is_bivariate (logical):  default := true
            %      selected_covariates (cell):  e.g., {'MergeAge', 'MergePtGender', 'MergeCdrsbBl', 'Components_1'}
            %      response (text):  one of selected_covariates, default := selected_covariates{end}.
            %      reuse_mdl (logical):  default true.
            
            ip = inputParser;
            addParameter(ip, "componentDir", pwd, @isfolder);
            addParameter(ip, "crossVal", 'on', @istext);
            addParameter(ip, "do_bayesopt", true, @islogical);
            addParameter(ip, "holdout", 0.2, @(x) isscalar && x < 1);
            addParameter(ip, "is_bivariate", true, @islogical);
            addParameter(ip, "selected_covariates", this.SELECTED_COVARIATES, @iscell);
            addParameter(ip, "response", "", @istext);
            addParameter(ip, "reuse_mdl", true, @islogical);
            parse(ip, varargin{:});
            ipr = ip.Results;
            this.componentDir = ipr.componentDir;
            assert(endsWith(this.componentDir, 'components'))
            this.crossVal = ipr.crossVal;
            this.do_bayesopt = ipr.do_bayesopt;
            this.holdout = ipr.holdout;
            this.is_bivariate = ipr.is_bivariate;
            this.selected_covariates = ipr.selected_covariates;
            this.response = ipr.response;
            if ~isempty(this.selected_covariates) && ...
                isemptytext(this.response)
                this.response = this.selected_covariates{end};
            end
            this.reuse_mdl = ipr.reuse_mdl;
        end
    end

    methods (Static)
        function t = beautify_varnames(t)
            % ADNIMerge compliant:
            names0 = {'AmyloidStatus' 'MergeAge' 'MergePtGender' 'MergePtEducat' 'MergePtEthCat' 'MergePtRacCat' 'MergeApoE4' 'MergeMmse' 'MergeDx' 'MergeCdrsbBl' 'CDMEMORY' 'CDORIENT' 'CDJUDGE' 'CDCOMMUN' 'CDHOME' 'CDCARE' 'CDGLOBAL' 'Phase' 'SITEID' 'ICV' 'PVE1'};
            % more readable in report:
            namesf = {'Amyloid Status' 'Age' 'Gender' 'Education' 'Ethnicity' 'Race' 'ApoE4' 'MMSE' 'Dx' 'CDR-SOB' 'CDR-Memory' 'CDR-Orientation' 'CDR-Judgement' 'CDR-Community' 'CDR-Home' 'CDR-Care' 'CDR-Global' 'ADNI Phase' 'Site ID' 'ICV' 'Gray'};
            
            matches = contains(t.Properties.VariableNames, names0);
            names0 = names0(matches);
            namesf = namesf(matches);
            t = renamevars(t, names0, namesf);
        end
        function t = combine_amyloid(t)
            %% center, scale, and average amyloid measures -> MergeAmyloid
            %  AmyloidStatus -> categories

            pib = normalize(t.MergePib); % pib centiloid -> z-score
            av45 = (196.9*t.MergeAv45) - 196.03; % see also "ADNI Centiloids Final.pdf"
            av45 = normalize(av45); % florbetapir -> centiloids -> z-score
            MergeAmyloid = mean([pib av45], 2, 'omitnan');
            t = addvars(t, MergeAmyloid, 'Before', 1);
            t = removevars(t, 'MergePib');
            t = removevars(t, 'MergeAv45');
            t = removevars(t, 'MergeAbeta');
            t = removevars(t, 'MergePibBl');
            t = removevars(t, 'MergeAv45Bl');
            t = removevars(t, 'MergeAbetaBl');

            arr = t.AmyloidStatus;
            arr(isnan(t.AmyloidStatus)) = 0;
            arr(1 == t.AmyloidStatus) = 1;
            arr(0 == t.AmyloidStatus) = -1;
            t.AmyloidStatus = arr;
        end
        function [u,v] = predictorImportance_fitrtree_pattern(indices, savefigs, close_fig, amy_status)
            %  Returns:
            %      u: table with ADNI variable names
            %      v: table with publishable variable names

            arguments
                indices {mustBeInteger,mustBeGreaterThan(indices,0)} = 1
                savefigs logical = false
                close_fig logical = false
                amy_status {mustBeNumeric} = []
            end

            if isempty(amy_status)
                ld = load('mladni_FDGDemographics_table_covariates.mat');
                t = ld.t;
                amy_stat = '';
                amy_stat_ = '';
            else
                if isnan(amy_status)
                    ld = load('mladni_FDGDemographics_table_covariates_amy_na.mat');
                    t = ld.amy_na;
                    amy_stat = 'Amyloid n/a';
                    amy_stat_ = 'amy_na';
                else
                    if logical(amy_status)
                        ld = load('mladni_FDGDemographics_table_covariates_amy_true.mat');
                        t = ld.amy_true;
                        amy_stat = 'Amyloid pos';
                        amy_stat_ = 'amy_true';
                    else
                        ld = load('mladni_FDGDemographics_table_covariates_amy_false.mat');
                        t = ld.amy_false;
                        amy_stat = 'Amyloid neg';
                        amy_stat_ = 'amy_false';
                    end
                end
            end

            % features, adding Component/pattern at far right of table
            u = mladni.NMFRegression.removenuisance(t);
            u = mladni.NMFRegression.removecomponents(u);
            u = u(:, [1 5 6 3 7   18 16 4 19 20   2 14 13 15 17   11 9 21 12 10   8   22 23 24 25 26 27 28]); % 28 vars remain
            u = mladni.NMFRegression.beautify_varnames(u);
 
            for idx = indices
                % reintroduce features, ending in Pattern := ld.t.Components(:,idx)
                Pattern = t.Components(:,idx); 
                v = addvars(u, Pattern); % added to far right of table
    
                Mdl = fitrtree(v, 'Pattern', 'PredictorSelection', 'curvature', 'Surrogate', 'on');
                imp = predictorImportance(Mdl);
                figure;
                bar(imp);
                title(sprintf('Fitted Regression Tree, %s Pattern %i', amy_stat, idx));
                ylabel('Predictor Importance Estimates', 'FontWeight', 'bold');
                xlabel('Predictors', 'FontWeight', 'bold');
                h = gca;
                h.XTickLabel = Mdl.PredictorNames;
                h.XTickLabelRotation = 45;
                h.TickLabelInterpreter = 'none';
                h.XTickLabelRotation = 90;
                len = length(Mdl.PredictorNames);
                set(gca, 'XTick', 1:len, 'XTickLabel', Mdl.PredictorNames)
            end

            if savefigs
                p = strcat('fitrtree_pattern_', amy_stat_);
                saveFigures(pwd, 'closeFigure', close_fig, 'prefix', p)
            end
        end
        function qp = queryPoint(t, pattern_index, cdr)
            arguments
                t table
                pattern_index {mustBeInteger,mustBePositive} = 5
                cdr double = 0.5
            end

            switch cdr
                case 0
                    t = t(t.CDGLOBAL == 0, :);
                case 0.5
                    t = t(t.CDGLOBAL == 0.5, :);
                case {1,2}
                    t = t(t.CDGLOBAL >= 0.5, :);
                otherwise
                    error('mladni:ValueError', 'NMFRegression.queryPoint.cdr->%g', cdr)
            end
            if ismember('Pattern', t.Properties.VariableNames)
                p = t.Pattern;
            end
            if ismember('Component', t.Properties.VariableNames)
                p = t.Component;
            end
            if ismember('Components', t.Properties.VariableNames)
                p = t.Components(:,pattern_index);
            end

            meanp = mean(p, 'omitnan');
            [~,idx] = min(abs(p - meanp));
            qp = t(idx,:);
        end
        function t = removecomponents(t)
            vns = t.Properties.VariableNames;
            hasc = contains(vns, 'Component');
            t = removevars(t, vns(hasc));
        end
        function t = removeids(t)
            t = removevars(t, {'Filelist' 'ImageDataID' 'MergeExamDate' 'EXAMDATE'});
        end
        function t = removenuisance(t)
            t = removevars(t, {'MergeDxBl' 'MergeTauBl' 'MergePTauBl' 'MergeFdgBl'}); % redundant
            t = mladni.NMFRegression.removeids(t);
            t = mladni.NMFRegression.siteid2categories(t);
            t = mladni.NMFRegression.combine_amyloid(t);
        end
        function t = siteid2categories(t)
            arr = t.SITEID;
            carr = num2cell(arr);
            carr = cellfun(@num2str, carr, 'UniformOutput', false);
            t.SITEID = carr;
        end
    end

    %% PROTECTED

    methods
        function [yFit,ySD,yInt] = kfoldPredict(this, varargin)
            %% web(fullfile(docroot, 'stats/regressiongam.predict.html'))
            %
            %  Params:
            %  Returns:
            %      vector of predicted response for predictor data in X.

            assert(strcmpi(this.crossVal, 'on'))

            ip = inputParser;
            ip.KeepUnmatched = true;
            parse(ip, varargin{:});
            ipr = ip.Results;

            cvmdl = this.mdl;
            try
                [yFit,ySD,yInt] = kfoldPredict(cvmdl, varargin{:});
            catch
                yFit = kfoldPredict(cvmdl, varargin{:});
            end
        end
        function mdl = mdl_bayesopt(this)
            if ~isempty(this.mdl_bayesopt_)
                mdl = this.mdl_bayesopt_;
                return
            end

            matfile = strcat(clientname(true, 2), '.mat');
            if isfile(matfile) && this.reuse_mdl
                ld = load(matfile);
                this.mdl_bayesopt_ = ld.mdl;
                mdl = this.mdl_bayesopt_;
                return
            end

            % web(fullfile(docroot, 'stats/fitrgam.html'))
            cvp = cvpartition(size(this.table_train, 1), 'KFold', 5);                
            maxNumSplits = optimizableVariable('maxNumSplits',[1,10],'Type','integer');
            numTrees = optimizableVariable('numTrees',[1,500],'Type','integer');
            minfun = @(z)kfoldLoss(fitrgam(this.table_train, this.response, ...
                'CVPartition', cvp, ...
                'MaxNumSplitsPerPredictor',z.maxNumSplits, ...
                'NumTreesPerPredictor',z.numTrees)); 
            results = bayesopt(minfun, [maxNumSplits, numTrees], ...
                'Verbose', 0, ...
                'IsObjectiveDeterministic', true, ...
                'AcquisitionFunctionName', 'expected-improvement-plus');
            zbest = bestPoint(results);
            mdl = fitrgam(this.table_train, this.response, ...
                'MaxNumSplitsPerPredictor', zbest.maxNumSplits, ...
                'NumTreesPerPredictor', zbest.numTrees, ...                
                'FitStandardDeviation', true, ...
                'OptimizeHyperparameters', 'all-univariate', ...
                'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                'ShowPlots', false, ...
                'Verbose', 0));
            save(matfile, 'mdl');
        end
        function mdl = mdl_cv(this)
            if ~isempty(this.mdl_cv_)
                mdl = this.mdl_cv_;
                return
            end

            matfile = strcat(clientname(true, 2), '.mat');
            if isfile(matfile) && this.reuse_mdl
                ld = load(matfile);
                this.mdl_cv_ = ld.mdl;
                mdl = this.mdl_cv_;
                return
            end

            % web(fullfile(docroot, 'stats/fitrgam.html'))
            mdl = fitrgam(this.table_train, ...
                this.response, ...
                'CrossVal', 'on');
            save(matfile, 'mdl');
        end
        function mdl = mdl_univariate(this)
            if ~isempty(this.mdl_uni_)
                mdl = this.mdl_uni_;
                return
            end

            matfile = strcat(clientname(true, 2), '.mat');
            if isfile(matfile) && this.reuse_mdl
                ld = load(matfile);
                this.mdl_uni_ = ld.mdl;
                mdl = this.mdl_uni_;
                return
            end

            mdl = fitrgam(this.table_train, ...
                this.response, ...
                'FitStandardDeviation', true, ...
                'OptimizeHyperparameters', 'all-univariate', ...
                'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                'ShowPlots', false, ...
                'Verbose', 0));
            save(matfile, 'mdl');
        end
        function mdl = mdl_bivariate(this)
            if ~isempty(this.mdl_bi_)
                mdl = this.mdl_bi_;
                return
            end            

            matfile = strcat(clientname(true, 2), '.mat');
            if isfile(matfile) && this.reuse_mdl
                ld = load(matfile);
                this.mdl_bi_ = ld.mdl;
                mdl = this.mdl_bi_;
                return
            end

            x = this.mdl_univariate.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
            mdl = fitrgam(this.table_train, ...
                this.response, ...
                'FitStandardDeviation', true, ...
                'InitialLearnRateForPredictors', x.InitialLearnRateForPredictors, ...
                'MaxNumSplitsPerPredictor', x.MaxNumSplitsPerPredictor, ...
                'NumTreesPerPredictor', x.NumTreesPerPredictor, ...
                'OptimizeHyperparameters', 'all-bivariate', ...
                'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName','expected-improvement-plus', ...
                'ShowPlots', false, ...
                'Verbose', 0));
            save(matfile, 'mdl');
        end
    end

    properties (Access = protected)
        cvp_
        compact_mdl_
        mdl_bayesopt_
        mdl_bi_
        mdl_cv_
        mdl_uni_
        table1pp_
        table_covariates_
        table_selected_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
