 classdef NMFRegression < handle
    %% line1
    %  line2
    %  
    %  Created 14-Sep-2022 20:29:21 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2039608 (R2022a) Update 5 for MACI64.  Copyright 2022 John J. Lee.
    
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
        component_weighted_average_csv
        idx_cn
        idx_dementia
        idx_mci
        table_covariates
        table_covariates_mat
        table_selected
        table_train
        table_test
    end

    methods

        %% GET

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
        function g = get.table_selected(this)
            %% last column is response variable

            if ~isempty(this.table_selected_)
                g = this.table_selected_;
                return
            end

            g = table;
            g_varnames = {};
            t = this.table_covariates;
            t_varnames = t.Properties.VariableNames;
            for v = t_varnames
                if contains(v{1}, this.selected_covariates)
                    g = addvars(g, t.(v{1}));
                    g_varnames = [g_varnames v{1}]; %#ok<AGROW> 
                end
            end
            g.Properties.VariableNames = g_varnames;
            g = splitvars(g, 'Components');
            try
                g = g(~matches(g.MergeDx, ''), :);
            catch
            end
            g = movevars(g, this.response, 'After', size(g,2));
        end
        function g = get.table_train(this)
            g = this.table_selected(training(this.cvpartition()),:);
        end
        function g = get.table_test(this)
            g = this.table_selected(test(this.cvpartition()),:);
        end

        %% 

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
        function mdl = fitrtree(this)
            mdl = fitrtree(this.table_train, this.response);
            view(mdl, "Mode", "graph")
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
        function this = regressionLearner(this)
            %% web(fullfile(docroot, 'stats/regressionlearner-app.html'))

            regressionLearner();  
        end

        function this = NMFRegression(varargin)
            %% NMFREGRESSION 
            %  Args:
            %      componentDir (folder):  containing component_weighted_average.csv
            %                              generatead by mladi.NMF.call(), 
            %                              e.g., baseline4/NumBases28/components
            %      crossVal (text):  default := 'on'.
            %      do_bayesopt (logical):  default := true
            %      holdout (scalar):  < 1 for training, testing with cross-validation.
            %      is_bivariate (logical):  default := true
            %      selected_covariates (cell):  e.g., {'Components', 'MergeAge', 'MergePtGender', 'MergeCdrsbBl'}
            %      response (text):  one of selected_covariates, default := selected_covariates{end}.
            %      reuse_mdl (logical):  default true.
            
            ip = inputParser;
            addParameter(ip, "componentDir", pwd, @isfolder);
            addParameter(ip, "crossVal", 'on', @istext);
            addParameter(ip, "do_bayesopt", true, @islogical);
            addParameter(ip, "holdout", 0.2, @(x) isscalar && x < 1);
            addParameter(ip, "is_bivariate", true, @islogical);
            addParameter(ip, "selected_covariates", {'Components', 'MergePtGender', 'MergeAge', 'MergeCdrsbBl'}, @iscell);
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
            if this.response == "" || isempty(this.response)
                this.response = this.selected_covariates{end};
            end
            this.reuse_mdl = ipr.reuse_mdl;
        end
    end

    %% PROTECTED

    properties (Access = protected)
        cvp_
        compact_mdl_
        mdl_bayesopt_
        mdl_bi_
        mdl_cv_
        mdl_uni_
        table_covariates_
        table_selected_
    end

    methods (Access = protected)
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
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
