classdef Test_NMFRegression < matlab.unittest.TestCase
    %% GAM Regrsssion with trees.
    %  
    %  Created 14-Sep-2022 20:29:22 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/test/+mladni_unittest.
    %  Developed on Matlab 9.12.0.2039608 (R2022a) Update 5 for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        componentDir
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mladni.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_fitrtree(this)
            this.testObj.fitrtree();
        end        
        function test_fsrftest(this)
            tbl = this.testObj.table_selected;
            head(tbl)
            
            [idx,scores] =  fsrftest(tbl, 'MergeCdrsbBl');
            vnames = strrep(tbl.Properties.VariableNames(idx),'_','\_');
            X = categorical(vnames);
            X = reordercats(X, vnames);
            bar(X, scores(idx))
            xlabel('Predictor rank')
            ylabel('Predictor importance score ~ -log(p)')
            %xticklabels(strrep(tbl.Properties.VariableNames(idx),'_','\_'))
            xtickangle(45)   

            idxInf = find(isinf(scores));
            if ~isempty(idxInf)
                fprintf('test_fsrftest:  idx of inf scores:');
                disp(idxInf)
                hold on
                bar(scores(idx(length(idxInf)+1))*ones(length(idxInf),1))
                legend('scores finite','scores < -eps')
                hold off
            end
        end
        function test_generalized_additive_model(this)
            t = this.testObj;
            response = t.table_train.MergeCdrsbBl;
            predictors = t.table_train.Components_3;
            shape_functions = [];
            this.testObj.generalized_additive_model( ...
                'response', response, ...
                'predictors', predictors, ...
                'shape_functions', shape_functions);
        end
        function test_histogramPartialDependence(this)
            this.testObj.histogramPartialDependence();
        end
        function test_mdls(this)
            %% web(fullfile(docroot, 'stats/train-generalized-additive-model-for-regression.html'))

            tic
            mdl = this.testObj.mdl(); % uses only tbl_training
            disp(mdl);
            CMdl = compact(mdl);
            whos('Mdl','CMdl')
            toc
            
            tbl_selected = this.testObj.table_selected;
            tbl_test = this.testObj.table_test;
            [yFit,ySD,yInt] = predict(CMdl,tbl_test);
            L = loss(CMdl,tbl_test);
            [yFit_nointeraction,ySD_nointeraction,yInt__nointeraction] = predict(CMdl,tbl_test,'IncludeInteractions',false);
            L_nointeractions = loss(CMdl,tbl_test,'IncludeInteractions',false);

            % Plot the sorted true responses together with the predicted responses and prediction intervals.
            yTrue = tbl_test.(this.testObj.response);
            [sortedYTrue,I] = sort(yTrue);
            
            figure 
            ax = nexttile;
            plot(sortedYTrue,'o')
            hold on
            plot(yFit(I))
            plot(yInt(I,1),'k:')
            plot(yInt(I,2),'k:')
            legend('True responses','Predicted responses', ...
                '95% Prediction interval limits','Location','best')
            title('Linear and interaction terms')
            hold off
            
            nexttile
            plot(sortedYTrue,'o')
            hold on
            plot(yFit_nointeraction(I))
            plot(yInt__nointeraction(I,1),'k:')
            plot(yInt__nointeraction(I,2),'k:')
            ylim(ax.YLim)
            title('Linear terms only')
            hold off

            % Interprete the prediction for the first test observation by using the plotLocalEffects function.
            % Predict a response value for the first observation of the test data, and plot the local effects of the terms in CMdl on the prediction. 
            % The predict function returns the prediction for the first observation table_test(1,:). 
            % The plotLocalEffects function creates a horizontal bar graph that shows the local effects of the terms in CMdl on the prediction. 
            % Each local effect value shows the contribution of each term to the prediction for table_test(1,:).
            yFit = predict(CMdl,tbl_test(1,:))

            figure
            numTerms = length(CMdl.PredictorNames);
            ple = plotLocalEffects(CMdl,tbl_test(1,:),'IncludeIntercept',true,'NumTerms',numTerms);
            ple.Parent.TickLabelInterpreter = 'none';

            % Compute the partial dependence values for MergeAge and plot the sorted values. 
            % Specify both the training and test data sets to compute the partial dependence values using both sets.
            [pd,x,y] = partialDependence(CMdl,'MergeAge', tbl_selected);
            [pd_sorted,I] = sort(pd);
            x_sorted = x(I);
            try
                x_sorted = reordercats(x_sorted,I);
            catch
            end
            figure
            plot(x_sorted,pd_sorted,'o:')
            xlabel('MergeAge')
            ylabel(this.testObj.response)
            title('Partial Dependence Plot')

            % Compute the partial dependence values for MergePtGender and plot the sorted values. 
            % Specify both the training and test data sets to compute the partial dependence values using both sets.
            [pd,x,y] = partialDependence(CMdl,'MergePtGender', tbl_selected);
            [pd_sorted,I] = sort(pd);
            x_sorted = x(I);
            try
                x_sorted = reordercats(x_sorted,I);
            catch
            end
            figure
            plot(x_sorted,pd_sorted,'o:')
            xlabel('MergePtGender')
            ylabel(this.testObj.response)
            title('Partial Dependence Plot')

        end
        function test_plot(this)
            plot(this.testObj);
        end
        function test_plotLocalEffects(this)
            ple = this.testObj.plotLocalEffects();
            disp(ple)
            x = ple.XEndPoints;
            y = ple.YEndPoints;
        end
        function test_plotPartialDependence(this)
            mdl = this.testObj.mdl();
            for p = mdl.PredictorNames
                figure
                plotPartialDependence(mdl, p{1}, "Conditional", "centered")
                ylabel('CDR sum-of-boxes baseline')
                xlabel(p{1})
            end
            saveFigures(pwd, 'closeFigure', false, 'prefix', 'individual conditional expectation Components_')
        end
        function test_predict(this)
            [yFit,ySD,yInt] = predict(this.testObj);
        end
        function test_regressionLearner(this)
            this.testObj_.regressionLearner();
            this.verifyClass(this.testObj.table_covariates, 'table');p,
        end
        function test_table_covariates(this)
            head(this.testObj.table_covariates)
        end
        function test_table_selected(this)
            disp(this.testObj.selected_covariates)
            head(this.testObj.table_selected)
        end
        function test_table_train(this)
            head(this.testObj.table_train)
            this.verifyEqual(size(this.testObj.table_train), [1245 31])
        end
        function test_table_test(this)
            head(this.testObj.table_test)
            this.verifyEqual(size(this.testObj.table_test), [138 31]);
        end
    end
    
    methods (TestClassSetup)
        function setupNMFRegression(this)
            import mladni.*
            this.componentDir = fullfile(getenv('SINGULARITY_HOME'), ...
                'ADNI', 'NMF_FDG', 'baseline4', 'NumBases32', 'components');
            this.testObj_ = NMFRegression('componentDir', this.componentDir);

            rng('default'); % for reproducibility
        end
    end
    
    methods (TestMethodSetup)
        function setupNMFRegressionTest(this)
            this.testObj = this.testObj_;
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
