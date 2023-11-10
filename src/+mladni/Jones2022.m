classdef Jones2022 < handle
    %% https://www.nature.com/articles/s41467-022-29047-4
    %
    %  "Relevant" objects are discussed in Jones' paper, esp. Table 2.
    %  
    %  Created 04-Mar-2023 11:44:49 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        N_PATTERNS = mladni.NMF.N_PATTERNS
    end

    properties (Dependent)
        componentsdir
        jones2022dir
        labels_check_mvregress % labels for plots
        labels_check_mvregress_2 % labels for plots, excluding NpBraak
        mask
        nbases
        plabels
        patt_path 
        relevant_file % mat file:  nii Filelist, merge/demographic variables, Components.  N = 781.
        responses_table2 % cell of var names
        responses_table2_2 % cell of var names, excluding NpBraak
        study_design
        workdir
    end

    methods % GET, SET
        function g = get.componentsdir(this)
            g = fullfile(this.workdir, 'baseline_cn', sprintf('NumBases%i', this.nbases), 'components');
        end
        function g = get.jones2022dir(~)
            g = fullfile(getenv('ADNI_HOME'), 'Jones_2022');
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
                sprintf('Jones2022_table_covariates_1stscan_%s.mat', this.study_design));
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
        function g_ = get.responses_table2_2(this)
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
        function g_ = get.labels_check_mvregress_2(this)
            g = this.labels_check_mvregress;
            g_ = g(~contains(g, 'Braak NFT'));
        end
        function g = get.study_design(this)
            g = this.study_design_;
        end
        function g = get.workdir(~)
            g = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
        end
    end

    methods
        function this = Jones2022(opts)
            arguments
                opts.study_design = 'longitudinal';
            end
            this.study_design_ = opts.study_design;
            if ~isfile(this.relevant_file)
                this.build_table_relevant();
            end
        end

        %% adjust spatial-autocorrelations

        function T = build_for_brainsmash(this)
            %% Prepares intermediates needed for brainsmash inference described in 
            %  test_neurodegeneration2.py/TestNeurodegeneration2.test_build_stats27.
            %  Requires availability of AFNI.3dmaskdump.
            %
            %  Returns table of 'SummaryTerm', 'TopicTermNumber', 'Filename'}

            workdir__ = fullfile(this.jones2022dir, 'Eigenbrains');
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
        function build_table_relevant(this)
            %% builds "Jones2022_table_covariates_1stscan_<study_design>.mat" using all visible baseline_* folders

            pwd0 = pushd(fileparts(this.relevant_file));
            
            ld = load(fullfile(this.componentsdir, ...
                sprintf('NMFCovariates_table_covariates_1stscan_%s.mat', this.study_design)));
            t = ld.t;

            % buld Components variable from Components_i
            vars = cellfun(@(x) sprintf('Components_%i', x), num2cell(1:this.N_PATTERNS), UniformOutput=false);
            mat = nan(size(t, 1), this.N_PATTERNS);
            for p = 1:this.N_PATTERNS
                mat(:,p) = t{:, vars{p}};
            end
            assert(~any(isnan(mat), 'all'))
            t = addvars(t, mat, NewVariableNames={'Components'});
            t = t(t.Cohort ~= categorical("CDR>0,amy-"), :);

            % save to filesystem
            save(this.relevant_file, 't');

            popd(pwd0);
        end
        function h = heatmap_eigenbrains(this, mat, clbls, rlbls, opts)
            %% selectd with FDR, but no considerations for spatial autocorrelations

            arguments
                this mladni.Jones2022
                mat double = []
                clbls cell = {}
                rlbls cell = {}
                opts.CellLabelFormat {mustBeTextScalar} = '%0.2g'
                opts.Colormap double = viridis
                opts.ColorScaling {mustBeTextScalar} = 'scaled' % 'scaledrows' 'scaledcolumns'
                opts.FlipColormap logical = false
                opts.FontSize {mustBeScalarOrEmpty} = 18
            end
            if isempty(mat)
                T = this.table();
                mat = T.corr;
                mat(T.fdr > 0.05) = 0;
            end
            if isempty(clbls)
                clbls = num2cell(1:this.nbases);
                clbls = cellfun(@(x) sprintf('P%i', x), clbls, UniformOutput=false);
            end
            if isempty(rlbls)
                rlbls = asrow(convertStringsToChars(T.term));
            end
            if opts.FlipColormap
                opts.Colormap = flipud(opts.Colormap);
            end

            sorted_bases_fn = fullfile(this.componentsdir, 'sorted_bases.mat');
            ld = load(sorted_bases_fn);
            sorted_bases = asrow(ld.sorted_bases);
            mat = mat(:, sorted_bases);

            figure;
            h = heatmap(clbls, rlbls, mat, ...
                CellLabelColor='auto', ...
                CellLabelFormat=opts.CellLabelFormat, ...
                Colormap=opts.Colormap, ColorScaling=opts.ColorScaling, ...
                FontSize=opts.FontSize);
            h.Title = "Pearson correlation of eigenbrains and ADNI patterns";
            h.YLabel = "Eigenbrains";
            h.XLabel = "ADNI Patterns";
        end
        function T = table_built_stats(this, varargin)
            %% test_neurodegeneration2.test_build_eigenbrains; 
            %  hand-assembled by importing test_build_eigenbrains.log

            if ~isempty(this.table_built_stats_)
                T = this.table_built_stats_;
                return
            end
            ld = load(fullfile(this.jones2022dir, 'testbuildeigenbrains.mat'));
            this.table_built_stats_ = ld.testbuildeigenbrains;
            T = this.table_built_stats_;            
            T = this.table_paren(T, varargin{:});
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

        %% table 2; mvregress related

        function h = heatmap_patterns_to_vars(this, sorted_bases_fn, dat_fn)
            arguments
                this mladni.Jones2022
                sorted_bases_fn {mustBeFile} = fullfile(this.componentsdir, 'sorted_bases.mat')
                dat_fn {mustBeFile} = fullfile(this.componentsdir, 'dat_20231108T160213.mat')
            end

            ld = load(sorted_bases_fn);
            sorted_bases = [1 asrow(ld.sorted_bases + 1)];
            ld = load(dat_fn);
            B = ld.dat.B';
            B = B(:, sorted_bases);
            clbls = ['intercept' this.plabels]; 
            rlbls = this.labels_check_mvregress_2;            

            h = this.heatmap(abs(B), clbls, rlbls, ColorScaling='scaledrows');
            title('abs(\beta)')
            ylabel('Dependent variables')
            xlabel('NMF pattern from ADNI CN')
            saveFigures(this.workdir, closeFigure=false, prefix=stackstr())
        end
        function dat = mvregress(this, maxiter)
            %% Following Jones' Table 2.
            %  Regression approach:
            %      - web(fullfile(docroot, 'stats/specify-the-response-and-study_design-matrices.html'))
            %      - d dimensions (Jones table-2 col-1) do share study_design matrix
            %      - when available, n observed subjects do share study_design matrix,
            %        but study_design may be incomplete for some subjects
            %  Multivariate GLM, Y = XB + E:
            %      - Y_{n x d}:  n subjects, d measured responses
            %      - X_{n x (p+1)}:  p predictors (NMF components), but sometime unavailable with up to
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
                maxiter double = 20000
            end

            ld = load(this.relevant_file);

            % responses from Jones' table 2
            Y = ld.t(:, this.responses_table2_2);
            Y = this.adjust_table_relevant(Y);
            sex = Y.Sex;
            sex = double(contains(sex, 'F'));
            sex(sex == 0) = -1;
            Y.Sex = sex;
            Y = Y{:,:};
            [n,d] = size(Y); % ~ 781 x 14

            X = zscore(ld.t.Components); % n x p ~ 781 x N_PATTERNS
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
            T1 = splitvars(T(:, 'Components')); % just the N_PATTERNS Components
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
                adjPval(r) = mdls{r}.ModelFitVsNullModel.Pvalue/this.N_PATTERNS;
                coeffs(r, :) = asrow(mdls{r}.Coefficients.Estimate);
            end

            % structured similarly to Jones' table 2
            T2 = table(ascol(this.labels_check_mvregress), N, R2, adjR2, preR2, AIC, coeffs(:,2:end), F, adjPval, ...
                VariableName={'Variable', 'N', 'R^2', 'adj R^2', 'pre R^2', 'AIC', '\beta P', 'F', 'adj p'});
            T2 = splitvars(T2);
        end

        %% bases compared

        function h = heatmap_patterns_to_eigenbrains(this)
            %% Kuehback-Leibler divergence, relative entropy from eigenbrains to patterns

            eb_path = '/Volumes/PrecunealSSD/Singularity/ADNI/Jones_2022/Eigenbrains';
            mat = nan(10, this.N_PATTERNS);

            for irow = 1:10

                eb = mlfourd.ImagingContext2(fullfile(eb_path, sprintf('MNI152_EB_%02i', irow)));    

                for icol = 1:this.N_PATTERNS        

                    patt = mlfourd.ImagingContext2(fullfile(this.patt_path, sprintf('Basis_%i.nii', icol))); 
                    div = eb.kldiv(patt, this.mask);    
                    mat(irow, icol) = div.nifti.img;
                end
            end

            sorted_bases_fn = fullfile(this.componentsdir, 'sorted_bases.mat');
            ld = load(sorted_bases_fn);
            sorted_bases = asrow(ld.sorted_bases);
            mat = mat(:, sorted_bases);

            clbls = this.plabels;
            rlbls = {'EB1' 'EB2' 'EB3' 'EB4' 'EB5' 'EB6' 'EB7' 'EB8' 'EB9' 'EB10'};
            h = mladni.Jones2022.heatmap(mat, clbls, rlbls, CellLabelFormat='%3.1f', ColorScaling='scaledrows', FlipColormap=false);
            ylabel('Eigenbrains');
            xlabel('ADNI Patterns');
            title('KL-Divergence of Patterns to Eigenbrains')
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
            colors = parula(length(cat));
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
                opts.Colormap double = viridis
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
                Colormap=opts.Colormap, ColorScaling=opts.ColorScaling, ...
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
        relevant_
        study_design_
        table_
        table_built_stats_
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
        function h = heatmap_beta(this, B)
            %B = B(2:end,:); % leave off intercept which can be large valued
            clbls = ['intercept' this.plabels];
            rlbls = this.labels_check_mvregress_2;
            figure;
            h = heatmap(B', XData=clbls, YData=rlbls, Colormap=parula, ColorScaling='scaledrows', FontSize=18);
            h.Colormap = parula;
            ylabel('Dependent Variables')
            xlabel('ADNI Pattern')
        end
        function h = heatmap_response(this, sigma)
            lbls = this.labels_check_mvregress_2;
            figure;
            h = heatmap(sigma, XData=lbls, YData=lbls, Colormap=parula, ColorScaling='log', FontSize=18);
            h.Colormap = parula;
            ylabel('Dependent Variables')
            xlabel('Dependent Variables')
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
