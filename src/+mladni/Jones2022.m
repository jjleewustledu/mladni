classdef Jones2022 < handle
    %% https://www.nature.com/articles/s41467-022-29047-4
    %
    %  "Relevant" objects are discussed in Jones' paper, esp. Table 2.
    %  
    %  Created 04-Mar-2023 11:44:49 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        N_PATTERNS = mladni.NMF.N_PATTERNS
        colormap = viridis
    end

    properties
        study_name
    end

    properties (Dependent)
        componentsdir
        data_home
        elabels
        labels_check_mvregress % labels for plots
        labels_check_mvregress_2 % labels for plots, excluding NpBraak
        mask
        nbases
        plabels
        patt_path 
        relevant_file % mat file:  nii Filelist, merge/demographic variables, Pattern_*.  N = 781.
        responses_table2 % cell of var names
        responses_table2_2 % cell of var names, excluding NpBraak
        sorted_bases
        study_design
        workdir
    end

    methods % GET, SET
        function g = get.componentsdir(this)
            g = fullfile(this.workdir, 'baseline_cn', sprintf('NumBases%i', this.nbases), 'components');
        end
        function g = get.data_home(~)
            g = fullfile(getenv('SINGULARITY_HOME'), 'MAYO');
        end
        function g = get.elabels(this)
            g = {'EB1' 'EB2' 'EB3' 'EB4' 'EB5' 'EB6' 'EB7' 'EB8' 'EB9' 'EB10'};
        end
        function g = get.mask(this)
            if ~isempty(this.mask_)
                g = this.mask_;
                return
            end
            mask_fqfn = fullfile(getenv('ADNI_HOME'), 'VolBin', 'mask.nii.gz');
            assert(isfile(mask_fqfn))
            this.mask_ = mlfourd.ImagingContext2(mask_fqfn);
            assert(all(size(this.mask_) == size(this.imaging{1})))
            g = this.mask_;
        end
        function     set.mask(this, s)
            this.mask_ = mlfourd.ImagingContext2(s);
        end        
        function g = get.nbases(this)
            g = this.N_PATTERNS;
        end
        function g = get.patt_path(this)
            g = fullfile(this.workdir, 'baseline_cn', sprintf('NumBases%i', this.nbases), 'OPNMF', 'niiImg');
        end
        function g = get.plabels(~)
            g = {'P1' 'P2' 'P3' 'P4' 'P5' 'P6' 'P7' 'P8' 'P9' 'P10' 'P11' 'P12'  ...
                 'P13' 'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21' 'P22' 'P23' 'P24'};
        end
        function g = get.relevant_file(this)
            g = fullfile(this.componentsdir, ...
                sprintf('NMFCovariates_table_covariates_1stscan_%s.mat', this.study_design));
        end
        function g = get.responses_table2(~)
            g = {'Metaroi', ...
                 'Age', 'MergeMmse', 'NpBraak', 'MergeCdrsb', ...
                 'META_TEMPORAL_SUVR', 'BRAAK1_SUVR', 'BRAAK34_SUVR', 'BRAAK56_SUVR', ...
                 'MergeHippocampus', 'CDGLOBAL', 'Sex', ...
                 'SUMMARYSUVR_WHOLECEREBNORM', ...
                 'ApoE4', 'ApoE2'};
            % BRAAK*:  tau PET
            % SUMMARSUVR*:  amyloid PET
        end
        function g_ = get.responses_table2_2(this)
            g = this.responses_table2;
            g_ = g(~contains(g, 'NpBraak', IgnoreCase=true));
        end
        function g = get.sorted_bases(this)
            if ~isempty(this.sorted_bases_)
                g = this.sorted_bases_;
                return
            end

            nmfh = mladni.NMFHierarchies();
            T = nmfh.table_patt_weighted_fdg();
            this.sorted_bases_ = asrow(T.indices_bases);
            g = this.sorted_bases_;
        end
        function g = get.labels_check_mvregress(~)
            g = {'FDG_{AD}', ...
                'Age', 'MMSE', 'Braak NFT', 'CDR-SOB', ...
                'Tau-PET' 'Tau-PET (Braak1)', 'Tau-PET (Braak34)', 'Tau-PET (Braak56)', ...
                'Hippo Vol', 'CDR', 'Sex', ...
                'Amyloid-PET', ...
                'E4+', 'E2+'};
        end
        function g_ = get.labels_check_mvregress_2(this)
            g = this.labels_check_mvregress;
            g_ = g(~contains(g, 'Braak NFT'));
        end
        function g = get.study_design(this)
            g = this.study_design_;
        end
        function g = get.workdir(this)
            g = fullfile(getenv('SINGULARITY_HOME'), this.study_name, 'NMF_FDG');
        end
    end

    methods
        function this = Jones2022(opts)
            arguments
                opts.study_design = 'longitudinal';
                opts.study_name = 'ADNI';  % 'MAYO'
            end
            this.study_design_ = opts.study_design;
            this.study_name = opts.study_name;
        end

        function T = build_for_brainsmash(this)
            %% Prepares intermediates needed for brainsmash inference described in 
            %  test_neurodegeneration2.py/TestNeurodegeneration2.test_build_stats27.
            %  Requires availability of AFNI.3dmaskdump.
            %
            %  Returns table of 'SummaryTerm', 'TopicTermNumber', 'Filename'}

            workdir__ = fullfile(this.data_home, 'Eigenbrains');
            cd(workdir__)

            Number = (1:10)';
            Label = cellfun(@(x) sprintf('EB%i', x), num2cell(Number), UniformOutput=false);
            Filename = cellfun(@(x) sprintf('MNI152_EB_%02i.nii.gz', x), num2cell(Number), UniformOutput=false);
            AltFileprefix = Label;
            T = table(Label, Number, Filename, AltFileprefix);
            parfor idx = 1:length(Number)
                system(sprintf('3dmaskdump -mask mask.nii.gz -o %s.txt -noijk %s', AltFileprefix{idx}, Filename{idx}))
            end
            writetable(T, sprintf('%s_T.csv', stackstr()));
        end
                
        function T = build_Pattern_var(this, T)
            assert(any(strcmp(this.table_relevant_.Properties.VariableNames, "Pattern_1")))
            assert(any(strcmp(this.table_relevant_.Properties.VariableNames, "Pattern_"+this.N_PATTERNS)))
            T = mergevars(T, ...
                "Pattern_" + (1:this.N_PATTERNS), ...
                NewVariableName="Pattern");
        end
        
        function T = table_built_stats(this, varargin)
            %% test_neurodegeneration2.test_build_eigenbrains; 
            %  hand-assembled by importing test_build_eigenbrains.log

            if ~isempty(this.table_built_stats_)
                T = this.table_built_stats_;
                return
            end
            ld = load(fullfile(this.data_home, 'neurodegeneration2_5k_EB.mat'));  % from brainsmash 5k
            this.table_built_stats_ = ld.T;
            T = this.table_built_stats_;            
            T = this.table_paren(T, varargin{:});
        end
        
        function t = table_relevant(this, varargin)
            if isempty(this.table_relevant_)
                ld = load(this.relevant_file);
                this.table_relevant_ = ld.t;
                if ~any(strcmp(this.table_relevant_.Properties.VariableNames, 'Pattern'))
                    this.table_relevant_ = this.build_Pattern_var(this.table_relevant_);
                end
            end
            t = this.table_relevant_;            
            t = mladni.AdniMerge.table_paren(t, varargin{:});
        end

        function T = table_termlist(this, varargin)
            %% enumerations of Jones Table 3

            if ~isempty(this.table_termlist_)
                T = this.table_termlist_;
                return
            end
            
            index = (1:10)';
            filename = cellfun(@(x) sprintf('EB%i.txt', x), num2cell(index), UniformOutput=false);            
            term = cellfun(@(x) sprintf('EB%i', x), num2cell(index), UniformOutput=false);
            this.table_termlist_ = table(index, filename, term);
            T = this.table_termlist_;            
            T = this.table_paren(T, varargin{:});
        end
        
        function T = table(this, varargin)
            %% compendium
            %  VariableNames:  'index', 'odd', 'term', 'corr', 'pval', 'fdr', 'groups'

            if ~isempty(this.table_)
                T = this.table_;
                return
            end

            % aufbau new table variables

            this.table_ = this.table_termlist;
            corr = zeros(length(this.table_.index), this.nbases);
            pval = ones(length(this.table_.index), this.nbases);
            fdr = ones(length(this.table_.index), this.nbases);
            this.table_ = addvars(this.table_, corr, pval, fdr, 'After', 'term', ...
                'NewVariableNames', {'corr', 'pval', 'fdr'});

            terms = unique(this.table_built_stats.term); % categorical
            for it = 1:length(terms)
                term_ = terms(it);

                % find idx of term_ in this.table_termlist
                it1 = find(this.table_.term == string(term_)); 
                if isempty(it1)
                    continue
                end
                assert(isscalar(it1))

                % term_-selected subtable of this.table_built_stats; 
                % avoid caches
                S = this.table_built_stats_(this.table_built_stats_.term == term_, ':'); 
                S = sortrows(S, 'basis');
                if this.nbases == length(S.basis)
                    r_row = S.r';
                    p_row = S.p';
                    this.table_.corr(it1,:) = r_row;
                    this.table_.pval(it1,:) = p_row;
                    continue
                end

                % special cases of U.basis
                r_row = nan(1, this.nbases);
                p_row = nan(1, this.nbases);
                for abasis = S.basis % int-valued
                    r_row(abasis) = S.r(abasis);
                    p_row(abasis) = S.p(abasis);
                end
                this.table_.corr(it1,:) = r_row;
                this.table_.pval(it1,:) = p_row;
            end
            
            % trim new table; update index -> index1
            %this.table_ = this.table_(~contains(this.table_.filename, ".txt"), ':');
            %this.table_ = this.table_(~any(isnan(this.table_.corr), 2), ':');
            %this.table_ = this.table_(~any(isnan(this.table_.pval), 2), ':');
            Nrows = size(this.table_, 1);
            this.table_ = addvars(this.table_, (1:Nrows)', 'After', 'index', 'NewVariableNames', 'index1');

            % do fdr
            for bi = 1:this.nbases
                [~,~,~,P] = fdr_bh(this.table_.pval(:,bi), 0.05, 'dep', 'yes');
                this.table_.fdr(:,bi) = ascol(P);
            end

            T = this.table_;
            T = this.table_paren(T, varargin{:});
        end

        %% heatmap of bases compared
        
        function h = heatmap_eigenbrains_pcorr(this)
            %% selectd with FDR, but no considerations for spatial autocorrelations
            
            % mat
            T = this.table();
            mat = T.corr;
            mat(T.fdr > 0.05) = 0;
            
            % clbls
            clbls = this.plabels;

            % rlbls
            rlbls = this.elabels;

            figure;
            this.heatmap(mat, clbls, rlbls);
            h.Title = "Pearson correlation, auto-correlations corrected";
            h.YLabel = "Eigenbrains";
            h.XLabel = "Patterns of covariance";
        end

        function h = heatmap_eigenbrains_kldiv(this)
            %% Kuehback-Leibler divergence, relative entropy from abs(eigenbrains) to patterns

            % mat
            eb_path = '/Volumes/PrecunealSSD/Singularity/MAYO/Eigenbrains';
            mat = nan(10, this.N_PATTERNS);
            for irow = 1:10
                eb = mlfourd.ImagingContext2(fullfile(eb_path, sprintf('MNI152_EB_%02i', irow)));   
                eb = abs(eb);
                for icol = 1:this.N_PATTERNS       
                    patt = mlfourd.ImagingContext2(fullfile(this.patt_path, sprintf('Basis_%i.nii', icol))); 
                    div = eb.kldiv(patt, this.mask);    
                    mat(irow, icol) = div.nifti.img;
                end
            end

            % clbls
            clbls = this.plabels;

            % rlbls
            rlbls = this.elabels;

            figure;
            h = this.heatmap(mat, clbls, rlbls, CellLabelFormat='%3.1f');
            ylabel('Eigenbrains');
            xlabel('Patterns of covariance');
            title('KL-Divergence')
        end

        %% table 2; mvregress related
        
        function dat = heatmaps_mvregress(this, maxiter)
            %% Following Jones' Table 2.
            %  Regression approach:
            %      - web(fullfile(docroot, 'stats/specify-the-response-and-study_design-matrices.html'))
            %      - d dimensions (Jones table-2 col-1) do share study_design matrix
            %      - when available, n observed subjects do share study_design matrix,
            %        but study_design may be incomplete for some subjects
            %  Multivariate GLM, Y = XB + E:
            %      - Y_{n x d}:  n subjects, d measured responses starting with Metaroi
            %      - X_{n x (p+1)}:  p predictors (NMF encodings), but sometime unavailable with up to
            %        p missing predictors, so use n-cell-array of {X^{(1..n)}}, 
            %        X^{(i)}_{d x K}:  K := d(p+1)
            %      - B_{(p+1) x d}
            %      - E_{n x d}
            %
            %  See also web(fullfile(docroot, 'stats/specify-the-response-and-study_design-matrices.html'))
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
                maxiter double = 100000
            end

            % responses from Jones' table 2
            T = this.table_relevant();
            Y = T(:, this.responses_table2_2);
            Y = this.adjust_table_relevant(Y);
            % sex = Y.Sex;
            % sex = double(contains(sex, 'F'));
            % Y.Sex = sex;
            Y = Y{:,:};
            [n,d] = size(Y); % ~ 781 x 14 vars from Jones' table 2

            X = zscore(this.table_relevant.Pattern); % n x p ~ 781 x N_PATTERNS
            p = size(X, 2);
            assert(size(X, 2) == p)
            Xmat = [ones(n,1) X];

            Xcell = cell(1,n);
            for i = 1:n
                Xcell{i} = [kron([Xmat(i,:)], eye(d))];
            end

            [beta,sigma,E,covB,logL] = mvregress(Xcell,Y, algorithm='ecm', varformat='beta', maxiter=maxiter);
            B = reshape(beta, d, p+1)';
            B = B(2:end, :);
            this.heatmap_beta(B)
            title('Multivariate regression \beta for dependent variables from encodings of patterns')            
            this.heatmap_response(sigma)
            title('Variance-covariance matrix of responses')
            this.set_gcf()

            %figure; pcolor(covB)
            %title('covB is variance-covariance of beta collapsed to vector')
            se = sqrt(diag(covB));
            SE = reshape(se, d, p+1)';
            SE = SE(2:end, :);
            this.heatmap_beta(SE)
            title('Standard errors')
            this.set_gcf()

            %[nlogL,COVB] = mvregresslike(Xcell, Y, beta, sigma, 'ecm');
            %figure; pcolor(COVB)
            %title('COVB')

            dat.beta = beta;
            dat.sigma = sigma;
            dat.E = E;
            dat.covB = covB;
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

            save( ...
                fullfile(this.componentsdir, sprintf('dat_%s.mat', datetime('now', 'Format','uuuuMMdd''T''HHmmss'))), ...
                'dat')
            saveFigures()
        end  

        function plot_relevant_all(this)
            this.plot_relevant('Metaroi', 'scatter', xlabel='FDG_{AD}');
            this.set_gcf()
            this.plot_relevant('MergeAge', 'scatter', xlabel='Age');
            this.set_gcf()
            this.plot_relevant('MergeMmse', 'scatter', xlabel='MMSE');
            this.set_gcf()
            this.plot_relevant('NpBraak', 'boxchart', catnames={'Stage 0' 'Stage 1' 'Stage 2' 'Stage 3' 'Stage 4' 'Stage 5' 'Stage 6'}, xlabel='Braak NFT', marksize=20);
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
            this.plot_relevant('CDGLOBAL', 'boxchart', catnames={'0' '0.5' '1' '2' '3'}, xlabel='CDR');
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
            c = T.Pattern;
            len = size(c, 2);

            h = figure; 
            tiledlayout(ceil(sqrt(len)), ceil(len/sqrt(len))); 
            try
                ylim_ = [min(c, [], "all"), max(c, [], "all")];
            catch ME
                fprintf("%s: %s has no ylim \n", stackstr(), tblvar);
                ylim_ = [];
            end
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
                if ~isempty(ylim_)
                    ylim(ylim_);
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
        
        %% table 2; fitlm related

        function [mdls,T2,h,h1] = fitlm(this)
            T = this.table_relevant();
            T1 = splitvars(T(:, 'Pattern')); % just the N_PATTERNS Components
            respname = this.responses_table2_2; % var names for T
            Nrows = length(respname);
            N = nan(Nrows, 1);
            R2 = nan(Nrows, 1);
            adjR2 = nan(Nrows, 1);
            % pre = nan(Nrows, 1);
            % sst = nan(Nrows, 1);
            % preR2 = nan(Nrows, 1);
            AIC = nan(Nrows, 1);
            BIC = nan(Nrows, 1);
            F = nan(Nrows, 1);
            adjPval = nan(Nrows, 1);            
            coeffs = nan(Nrows, size(T.Pattern, 2)+1);
            Ncoeffs = size(coeffs, 2);
            mdls = cell(size(respname));

            for r = 1:Nrows
                try
                    the_resp = T{:, respname{r}};
                    [the_resp,selected] = mladni.Jones2022.adjust_vars_relevant(the_resp, respname{r});
                    U1 = addvars(T1(selected,:), the_resp); % e.g., [Comp1, ..., Comp1, Metaroi]
                    mdls{r} = fitlm(U1);
                    N(r) = mdls{r}.NumObservations;
                    R2(r) = mdls{r}.Rsquared.Ordinary;
                    adjR2(r) = mdls{r}.Rsquared.Adjusted;
                    %arr = table2array(U1);
                    %arr = arr(~isnan(arr(:,end)),:);
                    %pre(r) = press(arr);
                    %sst(r) = sum(anova(mdls{r},'component',2).SumSq); % https://rpubs.com/RatherBit/102428#:~:text=adjusted%20R%2Dsquared%20%3D%201%20%2D,the%20lm%20and%20related%20functions.&text=So%2C%20you%20have%20to%20calculate,derive%20the%20predictive%20R%2Dsquared.
                    %preR2(r) = 1 - pre(r)/sst(r);
                    AIC(r) = mdls{r}.ModelCriterion.AIC;
                    BIC(r) = mdls{r}.ModelCriterion.BIC;
                    %F(r) = mdls{r}.ModelFitVsNullModel.Fstat;
                    adjPval(r) = mdls{r}.ModelFitVsNullModel.Pvalue/this.N_PATTERNS;
                    coeffs(r, :) = asrow(mdls{r}.Coefficients.Estimate(1:Ncoeffs));
                catch ME
                    handwarning(ME)
                end
            end

            % structured similarly to Jones' table 2
            T2 = table(ascol(this.labels_check_mvregress_2), N, R2, adjR2, AIC, BIC, coeffs(:,2:end), adjPval, ...
                VariableName={'Variable', 'N', 'R^2', 'adj R^2', 'AIC', 'BIC', '\beta P', 'adj p'}); %  preR2, 'pre R^2'
            T2 = splitvars(T2);
            disp(T2(:, [1:6, end]))

            beta = table2array(T2(:, 7:(end-1)))';
            abs_beta = abs(beta);
            h = this.heatmap_beta(beta, rlabels=this.labels_check_mvregress_2);
            title('Linear regression \beta for dependent variables from encodings of patterns')
            h1 = this.heatmap_beta(abs_beta, rlabels=this.labels_check_mvregress_2);
            title('Regression abs(\beta) for dependent variables from encodings of patterns')
        end
    end

    methods (Static)
        function T = adjust_table_relevant(T)
            for var = asrow(string(T.Properties.VariableNames))
                try
                    if contains(var, "Sex") 
                        T.(var) = double(strcmp(T.(var), 'M'));
                    else
                        v = T.(var);
                        T.(var) = (v - mean(v, 'omitmissing'))/std(v, 'omitmissing');
                    end
                    % range = max(T.(var)) - min(T.(var));
                    % assert(range > eps)
                    % T.(var) = (T.(var) - min(T.(var))) / range;
                catch ME
                    disp(ME.message)
                    fprintf("%s: rescaling is not appropriate for %s\n", stackstr(), var)
                end
            end
        end

        function [var,selected] = adjust_vars_relevant(var, respname)
            arguments
                var
                respname {mustBeTextScalar}
            end

            selected = true(size(var));

            if all(iscategorical(var))
                return
            end

            if contains(respname, "Sex")
                var = categorical(var);
                return
            end

            if contains(respname, "ApoE2")
                selected = ~isnan(var) & var >= 0;
                var = normalize(var(selected), "scale");
                return
            end

            if contains(respname, "Component")
                selected = ~isnan(var) & var >= 0;
                var = normalize(var(selected), "scale");
                return
            end

            selected = ~isnan(var) & var >= 0;
            var = normalize(var(selected), "scale");
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
        
        function [h,mu,sigma,q,notch] = al_goodplot(cat, catnames, metric)
            %% https://www.mathworks.com/matlabcentral/fileexchange/91790-al_goodplot-boxblot-violin-plot?s_tid=srchtitle_site_search_3_violin
            %  cat {mustBeNumeric}:  e.g., [2 4 6 8 ...]
            %  catnames {mustBeText}:  e.g., {'2' '4' '6' '8' ...}
            %  metric double:  cols have categories, rows have repetitions
            arguments
                cat {mustBeNumeric}
                catnames {mustBeText}
                metric double
            end

            % scrub cat if needed
            try
                if any(isnan(cat))
                    cat = cat(~isnan(cat));
                    metric = metric(:, ~isnan(cat));
                    assert(all(numel(cat) == size(metric,2))) % cols have categories, rows have repetitions
                end
            catch ME
                handwarning(ME)
            end

            % al_goodplot; see example.m in downloaded zip or Mathworks File Exchange
            % https://www.mathworks.com/matlabcentral/fileexchange/91790-al_goodplot-boxblot-violin-plot?s_tid=srchtitle_site_search_3_violin
            h = cell(size(cat));
            mu = nan(size(cat));
            sigma = nan(size(cat));
            q = cell(size(cat));
            notch = cell(size(cat));
            colors = this.colormap(length(cat));
            parzen = std(metric,0,"all")/1000;
            for cati = 1:length(cat)
                [h{cati},mu(cati),sigma(cati),q{cati},notch{cati}] = ...
                    al_goodplot(metric(:,cati), cat(cati), 0.618, colors(cati,:), 'bilateral', [], parzen, 0.309); 
            end
            xticks(cat);
            xticklabels(asrow(catnames));
            set(findobj(get(h{1}{1}, 'parent'), 'type', 'text'), 'fontsize', 18);
        end
        
        function h = heatmap(mat, clbls, rlbls, opts)
            arguments
                mat double
                clbls cell
                rlbls cell
                opts.CellLabelFormat {mustBeTextScalar} = '%0.2g'
                opts.Colormap double = cividis
                opts.ColorScaling {mustBeTextScalar} = 'scaled' % 'scaledrows' 'scaledcolumns'
                opts.FlipColormap logical = false
                opts.FontSize {mustBeScalarOrEmpty} = 18
            end
            if opts.FlipColormap
                opts.Colormap = flipud(opts.Colormap);
            end
            figure;
            h = heatmap(mat, ...
                XData=clbls, YData=rlbls, ...
                CellLabelFormat=opts.CellLabelFormat, ...
                Colormap=opts.Colormap, ...
                ColorScaling=opts.ColorScaling, ...
                FontSize=opts.FontSize);
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
        
        function t = table_paren(varargin)
            t = mladni.AdniMerge.table_paren(varargin{:});
        end
    end

    %% PROTECTED

    properties (Access = protected)
        mask_
        sorted_bases_
        study_design_
        table_
        table_built_stats_
        table_relevant_
        table_termlist_
    end

    methods (Access = protected)
        function hh = check_mvregress(this, dat)
            %  See also web(fullfile(docroot, 'stats/multivariate-general-linear-model.html'))

            lbls = this.labels_check_mvregress_2;            
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
        function h = heatmap_beta(this, B, opts)
            %B = B(2:end,:); % leave off intercept which can be large valued

            arguments
                this mladni.Jones2022
                B double
                opts.rlabels = this.labels_check_mvregress_2
            end

            clbls = this.plabels;
            rlbls = opts.rlabels;
            figure;
            h = heatmap(B', ...
                XData=clbls, YData=rlbls, ...
                CellLabelFormat='%0.2g', ...
                Colormap=this.colormap, ...
                ColorScaling='scaled', ...
                FontSize=18);
            h.Colormap = this.colormap;
            ylabel('Dependent variables')
            xlabel('Encodings of patterns')
        end
        function h = heatmap_response(this, sigma, opts)
            arguments
                this mladni.Jones2022
                sigma double
                opts.labels = this.labels_check_mvregress_2
            end

            lbls = opts.labels;
            figure;
            h = heatmap(sigma, ...
                XData=lbls, YData=lbls, ...
                CellLabelFormat='%0.2g', ...
                Colormap=this.colormap, ...
                ColorScaling='log', ...
                FontSize=18);
            h.Colormap = this.colormap;
            ylabel('Dependent Variables')
            xlabel('Dependent Variables')
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
