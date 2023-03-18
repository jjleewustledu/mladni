classdef Jones2022
    %% https://www.nature.com/articles/s41467-022-29047-4
    %
    %  "Relevant" objects are discussed in Jones' paper, esp. Table 2.
    %  
    %  Created 04-Mar-2023 11:44:49 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Dependent)
        labels_check_mvregress % labels for plots
        labels_check_mvregress_ % labels for plots, excluding NpBraak
        relevant_file % mat file:  nii Filelist, merge/demographic variables, Components
        responses_table2 % cell of var names
        responses_table2_ % cell of var names, excluding NpBraak
    end

    methods % GET
        function g = get.relevant_file(~)
            g = fullfile(getenv('ADNI_HOME'), 'NMF_FDG', ...
                'mladni_FDGDemographics_table_covariates_on_cn_relevant_baselines.mat');
        end
        function g = get.responses_table2(~)
            g = {'Metaroi', ...
                 'MergeAge', 'MergeMmse', 'NpBraak', 'MergeCdrsb', ...
                 'META_TEMPORAL_SUVR', 'BRAAK1_SUVR', 'BRAAK34_SUVR', 'BRAAK56_SUVR', ...
                 'MergeHippocampus', 'CDGLOBAL', 'Sex', ...
                 'SUMMARYSUVR_WHOLECEREBNORM', ...
                 'MergeApoE4', 'ApoE2'};
            % BRAAK*:  tau PET
            % SUMMARSUVR*:  amyloid PET
        end
        function g_ = get.responses_table2_(this)
            g = this.responses_table2;
            g_ = g(~contains(g, 'NpBraak', IgnoreCase=true));
        end
        function g = get.labels_check_mvregress(~)
            g = {'FDG_{AD}', ...
                'Age', 'MMSE', 'Braak NFT', 'CDR-SOB', ...
                'Tau-PET' 'Tau-PET (Braak1)', 'Tau-PET (Braak34)', 'Tau-PET (Braak56)', ...
                'Hippo Vol', 'CDR', 'Sex', ...
                'Amyloid-PET', ...
                'E4+', 'E2+'};
        end
        function g_ = get.labels_check_mvregress_(this)
            g = this.labels_check_mvregress;
            g_ = g(~contains(g, 'Braak NFT'));
        end
    end

    methods
        function this = Jones2022(varargin)
        end
        function build_table_relevant(this)
            pwd0 = pushd(fileparts(this.relevant_file));
            
            g = glob('baseline_*');
            g = g(~contains(g, 'cdr_ge'));
            ld = load(fullfile(g{1}, 'NumBases16', 'components', 'mladni_FDGDemographics_table_covariates_on_cn.mat'));
            t = ld.t;
            for i = 2:length(g)
                ld = load(fullfile(g{i}, 'NumBases16', 'components', 'mladni_FDGDemographics_table_covariates_on_cn.mat'));
                t = [t; ld.t];
            end
            save(this.relevant_file, 't');

            popd(pwd0);
        end

        %% table 2; mvregress related

        function dat = mvregress(this, maxiter)
            %% Following Jones' Table 2.
            %  Regression approach:
            %      - web(fullfile(docroot, 'stats/specify-the-response-and-design-matrices.html'))
            %      - d dimensions (Jones table-2 col-1) do share design matrix
            %      - when available, n observed subjects do share design matrix,
            %        but design may be incomplete for some subjects
            %  Multivariate GLM, Y = XB + E:
            %      - Y_{n x d}:  n subjects, d measured responses
            %      - X_{n x (p+1)}:  p predictors (NMF components), but sometime unavailable with up to
            %        p missing predictors, so use n-cell-array of {X^{(1..n)}}, 
            %        X^{(i)}_{d x K}:  K := d(p+1)
            %      - B_{(p+1) x d}
            %      - E_{n x d}
            %
            %  See also web(fullfile(docroot, 'stats/specify-the-response-and-design-matrices.html'))
            %  See also web(fullfile(docroot, 'stats/multivariate-general-linear-model.html'))

            % dat = j.mvregress(8*800)
            % Warning: Maximum iterations completed. Convergence criterion not satisfied. 
            % > In mvregress (line 494)
            % In mladni/Jones2022/mvregress (line 164)   

            % dat = j.mvregress(51200)
            % Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  2.571809e-18. 
            % > In mvregress (line 525)
            % In mladni/Jones2022/mvregress (line 164)

            arguments
                this mladni.Jones2022
                maxiter double = 20000
            end

            ld = load(this.relevant_file);

            % responses from Jones' table 2
            Y = ld.t(:, this.responses_table2_);
            Y = this.adjust_table_relevant(Y);
            sex = Y.Sex;
            sex = double(contains(sex, 'F'));
            sex(sex == 0) = -1;
            Y.Sex = sex;
            Y = Y{:,:};
            [n,d] = size(Y); % ~ 781 x 14

            X = zscore(ld.t.Components); % n x p ~ 781 x 16
            p = size(X, 2);
            assert(size(X, 2) == p)
            Xmat = [ones(n,1) X];

            Xcell = cell(1,n);
            for i = 1:n
                Xcell{i} = [kron([Xmat(i,:)], eye(d))];
            end

            [beta,sigma,E,V,logL] = mvregress(Xcell,Y, algorithm='ecm', varformat='beta', maxiter=maxiter);
            B = reshape(beta, d, p+1)';
            this.heatmap_beta(B) % excludes large intercepts for age, mmse
            title('\beta matrix')
            this.heatmap_response(sigma)
            title('variance-covariance matrix of responses')
            this.set_gcf()

            %figure; pcolor(V)
            %title('V is variance-covariance of beta collapsed to vector')
            se = sqrt(diag(V));
            SE = reshape(se, d, p+1)';
            this.heatmap_beta(SE)
            title('standard errors')
            this.set_gcf()

            %[nlogL,COVB] = mvregresslike(Xcell, Y, beta, sigma, 'ecm');
            %figure; pcolor(COVB)
            %title('COVB')

            dat.beta = beta;
            dat.sigma = sigma;
            dat.E = E;
            dat.V = V;
            dat.logL = logL;
            dat.se = se;
            dat.B = B;
            dat.n = n;
            dat.d = d;
            dat.p = p;
            %dat.nlogL = nlogL;
            %dat.COVB = COVB;

            this.check_mvregress(dat)
            this.set_gcf()

            save(sprintf('dat_%s.mat', datetime('now', 'Format','uuuuMMdd''T''HHmmss')), 'dat')
        end     
        function h = plot_results(this, dat)
            %% plot dat after mvregress

            h = figure;
        end

        function plotall_relevant(this)
            this.plot_relevant('Metaroi', 'scatter', xlabel='FDG_{AD}');
            this.set_gcf()
            this.plot_relevant('MergeAge', 'scatter', xlabel='Age');
            this.set_gcf()
            this.plot_relevant('MergeMmse', 'scatter', xlabel='MMSE');
            this.set_gcf()
            this.plot_relevant('NpBraak', 'boxchart', catnames={'Stage 0' 'Stage 1' 'Stage 2' 'Stage 3' 'Stage 5' 'Stage 6'}, xlabel='Braak NFT', marksize=20);
            this.set_gcf()
            this.plot_relevant('MergeCdrsb', 'scatter', xlabel='CDR-SOB');
            this.set_gcf()
            this.plot_relevant('META_TEMPORAL_SUVR', 'scatter', xlabel='Tau-PET');
            this.set_gcf()
            this.plot_relevant('BRAAK1_SUVR', 'scatter', xlabel='Tau-PET (Braak1)');
            this.set_gcf()
            this.plot_relevant('BRAAK34_SUVR', 'scatter', xlabel='Tau-PET (Braak34)');
            this.set_gcf()
            this.plot_relevant('BRAAK56_SUVR', 'scatter', xlabel='Tau-PET (Braak56)');
            this.set_gcf()
            this.plot_relevant('MergeHippocampus', 'scatter', xlabel='Hippo Vol');
            this.set_gcf()
            this.plot_relevant('CDGLOBAL', 'boxchart', catnames={'0' '1' '1.5' '2'}, xlabel='CDR');
            this.set_gcf()
            this.plot_relevant('Sex', 'boxchart', catnames={'Female','Male'}, xlabel='Sex');
            this.set_gcf()
            this.plot_relevant('SUMMARYSUVR_WHOLECEREBNORM', 'scatter', xlabel='Amyloid-PET');
            this.set_gcf()
            this.plot_relevant('MergeApoE4', 'boxchart', catnames={'0 alleles' '1 allele' '2 alleles'}, xlabel='E4+');
            this.set_gcf()
            this.plot_relevant('ApoE2', 'boxchart', catnames={'0 alleles' '1 allele' '2 alleles'}, xlabel='E2+');
            this.set_gcf()
        end
        function h = plot_relevant(this, tblvar, plotfun, opts)
            %% Following Jones' Table 2, mapping tblvar -> averaged components
            %  Args:
            %     this mladni.Jones2022
            %     tblvar {mustBeTextScalar} = 'Metaroi' 
            %     'MergeAge', 'MergeMmse', 'NpBraak', 'MergeCdrsb', 
            %     'BRAAK1_SUVR', 'BRAAK34_SUVR', 'BRAAK56_SUVR' (tau-pet), 'MergeHippocampus', 'CDGLOBAL', 'Sex', 
            %     'SUMMARYSUVR_WHOLECEREBNORM' (amyloid-pet), 'MergeApoE4', 'ApoE2'
            %     plotfun {mustBeTextScalar} = 'scatter'
            %     opts.tbl table = this.table_relevant()
            %     opts.catnames cell = {}
            %     opts.xlabel {mustBeTextScalar} = "FDG_{AD}"
            %     opts.ylabel {mustBeTextScalar} = "\beta_"
            %     opts.marksize double = 4 % boxcharts only

            arguments
                this mladni.Jones2022
                tblvar {mustBeTextScalar} = 'Metaroi'
                plotfun {mustBeTextScalar} = 'scatter'
                opts.tbl table = this.table_relevant()
                opts.catnames cell = {}
                opts.xlabel {mustBeTextScalar} = "FDG_{AD}"
                opts.ylabel {mustBeTextScalar} = "ADNI"
                opts.marksize double = 4 % boxcharts only
            end
            if isempty(opts.catnames)
                opts.catnames = opts.tbl.(tblvar);
            end

            T = opts.tbl;
            c = T.Components;
            len = size(c, 2);

            h = figure; 
            tiledlayout(upper(sqrt(len)), upper(len/sqrt(len))); 
            for ci = 1:len
                nexttile; 
                if strcmp(plotfun, "boxchart")
                    hh = this.boxchart(T.(tblvar), opts.catnames, c(:,ci), opts.marksize);
                else
                    try
                        hh = feval(plotfun, T.(tblvar), c(:,ci), ".");
                    catch
                        hh = this.boxchart(T.(tblvar), opts.catnames, c(:,ci), opts.marksize);                        
                    end
                end
                set(gca, 'FontSize', 16);

                %title(sprintf("ADNI P%i", ci));
                if contains(opts.ylabel, "ADNI")
                    the_ylabel = opts.ylabel+" P"+ci;
                end
                ylabel(the_ylabel, FontSize=20, FontWeight="bold");
                xlabel(opts.xlabel, FontSize=20, FontWeight="bold");
            end
        end
        function t = table_relevant(this, varargin)
            if isempty(this.relevant_)
                ld = load(this.relevant_file);
                this.relevant_ = ld.t;
            end
            t = this.relevant_;            
            t = mladni.AdniMerge.table_paren(t);
        end
        
        %% table 2; fitlm related

        function [mdls,T2] = fitlm(this)
            T = this.table_relevant();
            T = this.adjust_table_relevant(T);
            T1 = splitvars(T(:, 'Components')); % just the 16 Components
            for col = 1:size(T1,2)
                T1{:,col} = zscore(T1{:,col});
            end
            respname = this.responses_table2; % var names for T
            Nrows = length(respname);
            N = nan(Nrows, 1);
            R2 = nan(Nrows, 1);
            adjR2 = nan(Nrows, 1);
            pre = nan(Nrows, 1);
            sst = nan(Nrows, 1);
            preR2 = nan(Nrows, 1);
            AIC = nan(Nrows, 1);
            BIC = nan(Nrows, 1);
            F = nan(Nrows, 1);
            adjPval = nan(Nrows, 1);            
            coeffs = nan(Nrows, size(T.Components, 2)+1);

            parfor r = 1:Nrows
                the_resp = T{:, respname{r}};
                if strcmp(respname{r}, 'Sex')
                    sex_ = double(strcmp(the_resp, 'F'));
                    sex_(sex_~=1) = -1;
                    the_resp = sex_;
                end
                %the_resp = zscore(the_resp);
                U1 = addvars(T1, the_resp); % e.g., [Comp1, ..., Comp1, Metaroi]
                mdls{r} = fitlm(U1); 
                N(r) = mdls{r}.NumObservations;
                R2(r) = mdls{r}.Rsquared.Ordinary;
                adjR2(r) = mdls{r}.Rsquared.Adjusted;
                arr = table2array(U1);
                arr = arr(~isnan(arr(:,end)),:);
                pre(r) = press(arr);
                sst(r) = sum(anova(mdls{r},'component',2).SumSq); % https://rpubs.com/RatherBit/102428#:~:text=adjusted%20R%2Dsquared%20%3D%201%20%2D,the%20lm%20and%20related%20functions.&text=So%2C%20you%20have%20to%20calculate,derive%20the%20predictive%20R%2Dsquared.
                preR2(r) = 1 - pre(r)/sst(r);
                AIC(r) = mdls{r}.ModelCriterion.AIC;
                BIC(r) = mdls{r}.ModelCriterion.BIC;
                F(r) = mdls{r}.ModelFitVsNullModel.Fstat;
                adjPval(r) = mdls{r}.ModelFitVsNullModel.Pvalue/16;
                coeffs(r, :) = asrow(mdls{r}.Coefficients.Estimate);
            end

            % structured similarly to Jones' table 2
            T2 = table(ascol(this.labels_check_mvregress), N, R2, adjR2, preR2, AIC, coeffs(:,2:end), F, adjPval, ...
                VariableName={'Variable', 'N', 'R^2', 'adj R^2', 'pre R^2', 'AIC', '\beta P', 'F', 'adj p'});
            T2 = splitvars(T2);
        end

    end

    methods (Static)
        function T = adjust_table_relevant(T)
            T.MergeAge = T.MergeAge/100; % centuries
            T.MergeMmse = T.MergeMmse/30; % fractional
            T.MergeHippocampus = T.MergeHippocampus/1e4;
        end        
        function h = boxchart(cats, catnames, comp, marksize)
            try
                if any(isnan(cats))
                    cats = cats(~isnan(cats));
                    comp = comp(~isnan(cats));
                end
            catch
            end
            x = categorical(cats, sort(unique(cats)), catnames);
            h = boxchart(x, comp, Notch="on");
            set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontsize', 18);
            hold on
            s = swarmchart(x, comp, marksize, MarkerFaceAlpha=0.5, MarkerEdgeAlpha=0.5);
            %s.XJitter = 'rand';
            %s.XJitterWidth = 0.25;
            hold off
        end
        function errors = regf(X1train,X2train,ytrain,X1test,X2test,ytest)
            tbltrain = table(X1train,X2train,ytrain, ...
                'VariableNames',{'Acceleration','Displacement','Weight'});
            tbltest = table(X1test,X2test,ytest, ...
                'VariableNames',{'Acceleration','Displacement','Weight'});
            mdl = fitlm(tbltrain,'Weight ~ Acceleration + Displacement');
            yfit = predict(mdl,tbltest);
            MAE = mean(abs(yfit-tbltest.Weight));
            adjMAE = MAE/range(tbltest.Weight);
            errors = [MAE adjMAE];
        end
        function set_gcf(opts)
            arguments
                opts.fracx double = 0.5
                opts.Npx double = 3400
                opts.fracy double = 1
                opts.Npy double = 1440
            end
            set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
        end
    end

    %% PROTECTED

    properties (Access = protected)
        relevant_
    end

    methods (Access = protected)
        function hh = check_mvregress(this, dat)
            %  See also web(fullfile(docroot, 'stats/multivariate-general-linear-model.html'))

            lbls = this.labels_check_mvregress_;            
            hh = figure();

            tiledlayout(3,5);
            for tile = 1:14
                nexttile;
                            
                z = dat.E/chol(dat.sigma);
                plot(z(:,1),z(:,tile),'.')
                xlabel(lbls{1}, FontSize=24, FontWeight="bold")
                ylabel(lbls{tile}, FontSize=24, FontWeight="bold")
                hold on

                % Overlay standard normal contours
                z1 = linspace(-5,5);
                z2 = linspace(-5,5);
                [zx,zy] = meshgrid(z1,z2);
                zgrid = [reshape(zx,100^2,1),reshape(zy,100^2,1)];
                zn = reshape(mvnpdf(zgrid),100,100);
                [c,h] = contour(zx,zy,zn);
                clabel(c,h, FontSize=12, FontWeight="normal")

                hold off
            end
        end   
        function h = heatmap_beta(this, B)
            %B = B(2:end,:); % leave off intercept which can be large valued
            clbls = {'intercept' 'P1' 'P2' 'P3' 'P4' 'P5' 'P6' 'P7' 'P8' 'P9' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16'};
            rlbls = this.labels_check_mvregress_;
            figure;
            h = heatmap(B', XData=clbls, YData=rlbls, Colormap=parula, ColorScaling='scaledrows', FontSize=18);
            h.Colormap = parula;
            ylabel('Dependent Variables')
            xlabel('ADNI Pattern')
        end
        function h = heatmap_response(this, sigma)
            lbls = this.labels_check_mvregress_;
            figure;
            h = heatmap(sigma, XData=lbls, YData=lbls, Colormap=parula, ColorScaling='log', FontSize=18);
            h.Colormap = parula;
            ylabel('Dependent Variables')
            xlabel('Dependent Variables')
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
