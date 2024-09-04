classdef NMFRadar < handle
    %% Warning: For the log scale values, recommended axes limit is [1.000000e-04, 10] with an axes interval of 5. 
    %  For managing natural ordering of labels, filename, etc., see also
    %  https://blogs.mathworks.com/pick/2014/12/05/natural-order-sorting/
    %  https://www.mathworks.com/matlabcentral/fileexchange/34464-customizable-natural-order-sort?s_tid=prof_contriblnk
    %  https://www.mathworks.com/matlabcentral/fileexchange/47433-natural-order-row-sort?s_tid=prof_contriblnk
    %  https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort?s_tid=prof_contriblnk
    %  
    %  Created 16-Feb-2023 22:57:16 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        AXES_COLOR = [0.9 0.9 0.9];
        N_PATTERNS = mladni.NMF.N_PATTERNS
        ERR_COLOR = [0.9 0.9 0.9]
        ERR_ALPHA = 0.08
        GROUP_COLORS = [
            0.134633333333333 0.240916666666667 0.425183333333333;
            0.338566666666667 0.361233333333333 0.4268;
            0.4879 0.4841 0.47105;
            0.654266666666667 0.616366666666667 0.4597;
            0.833733333333333 0.7609 0.393416666666667]  % cividis(5)
        GROUP_LINEWIDTHS = [1,1,1,1,1]
    end

    properties
        groups = { ...
            '' 'CDR=0,amy+' 'CDR=0.5,amy+' 'CDR>0.5,amy+'} 
        groupLabels = {...
            '' 'CDR=0,amy+' 'CDR=0.5,amy+' 'CDR>0.5,amy+'}
        groups0 = { ...
            'cn', ...
            'preclinical', ...
            'cdr_0p5_apos', ...
            'cdr_gt_0p5_apos'}
        mergeDx = { ...
            'CDR=0,amy-' 'CDR=0,amy+' 'CDR=0.5,amy+' 'CDR>0.5,amy+'}
        
        matfile0 = 'NMFCovariates_table_cn_1stscan_longitudinal.mat'
        matfile_cohort = 'patterns_of_neurodegeneration_20240630_for_import.mat' 
        % Output from R: Singularity/ADNI/NMF_FDG/patterns_of_neurodegeneration_20230921.Rmd
        % b2 <- gam(list(
        % y0~s(age,k=20,by=interaction(sex))+sex+apoe4+cohort, y1~s(age,k=20,by=interaction(sex))+sex+apoe4+cohort, y2~s(age,k=20,by=interaction(sex))+sex+apoe4+cohort, ...
        % family=mvn(d=24), data=soto)
        % b2
        % summary(b2)  
        matfile_cohort_sorted = 'CohortCoefficientsSorted.mat'  % P1 has highest SUVR, P24 has lowest SUVR
        workdir
    end

    properties (Dependent)
        AxesLabels
        AxesLabelsNull
        figdir
        N_bases_target
        N_groups
        sorted_bases % indexing from VolBin -> indexing sorted by pattern-weighted FDG SUvR; 
                     % e.g., sorted_bases(16) == 20
    end

    methods % GET
        function g = get.AxesLabelsNull(~)
            g = {'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''};
        end
        function g = get.AxesLabels(this)
            g = {'p_1', 'p_2', 'p_3', 'p_4', ...
                 'p_5', 'p_6', 'p_7', 'p_8', ...
                 'p_9', 'p_{10}', 'p_{11}', 'p_{12}', ...
                 'p_{13}', 'p_{14}', 'p_{15}', 'p_{16}', ...
                 'p_{17}', 'p_{18}', 'p_{19}', 'p_{20}', ...
                 'p_{21}', 'p_{22}', 'p_{23}', 'p_{24}'};
            g = g(this.sorted_bases);
        end
        function g = get.sorted_bases(this)
            if ~isempty(this.sorted_bases_)
                g = this.sorted_bases_;
                return
            end
            T = this.table_patt_weighted_fdg();
            g = asrow(T.indices_bases);
        end
        function g = get.figdir(this)
            g = fullfile(this.workdir, 'baseline_cn', 'results');
        end        
        function g = get.N_bases_target(this)
            g = this.N_PATTERNS;
        end
        function g = get.N_groups(this)
            g = length(this.groups);
        end
    end

    methods
        function [s1,s2] = call_intercept(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
                opts.pvalues logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            PC = ld.patternsofneurodegeneration20240630forimport;
            intercept = PC(contains(PC.ParamCoefficients, "Intercept"), :);
            if opts.show
                disp(intercept)
            end

            % Estimate
            P = intercept.Estimate';
            CI = 1.96*intercept.StdErr';
            amin = min(P-CI, [], 'all');
            amax = max(P+CI, [], 'all');
            h = figure;
            c = this.suvr2viridis(min(P, [], "all"), suvr_min=0.76, suvr_max=1.36);
            s1 = plot_with_CI(this, P, CI, AxesMin=amin, AxesMax=amax, Color=c);
            saveFigure2(h, fullfile(this.figdir, 'call_intercept_GAM_beta_intercept'));
    
            if opts.pvalues
                % PValue
                fprintf('%s:\n', stackstr())
                [h, crit_p, adj_ci_cvrg, Padj] = fdr_bh(intercept.PValue', 0.05, 'dep', 'yes');
                fprintf("crit_p->%g; adj_ci_cvrg->%g", crit_p, adj_ci_cvrg)
                %disp(this.AxesLabels(h))
                disp(ascol(Padj(this.sorted_bases)))
                Padj = asrow(Padj); 
                amin = min(Padj, [], 'all');
                amax = max(Padj, [], 'all');
                if amin < 1e-15 || amax < 1e-15
                    return
                end
                figure
                axes_scaling = repmat({'log'}, [1 this.N_PATTERNS]);
                s2 = plot(this, Padj, AxesMin=amin, AxesMax=amax, ...
                    legend={'CDR=0, amy-'}, ti='FDR p-value \beta_{Intercept}', AxesScaling=axes_scaling);
                % saveFigures(this.figdir, closeFigure=false, prefix='call_intercept_FDR_p-value_beta_intercept');
            end
        end

        function [s1,s2] = call_apoe4(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
                opts.pvalues logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            PC = ld.patternsofneurodegeneration20240630forimport;
            apoe4 = PC(contains(PC.ParamCoefficients, "apoe4"), :);
            if opts.show
                disp(apoe4)
            end

            % Estimate
            %load("P.mat");
            %load("SE.mat");
            P = apoe4.Estimate';
            CI = 1.96*apoe4.StdErr';
            amin = min(P-CI, [], 'all');
            amax = max(P+CI, [], 'all');
            figure
            c = this.suvr2viridis(min(P, [], "all"), suvr_min=-0.254, suvr_max=0.0108);
            s1 = plot_with_CI(this, P, CI, AxesMin=amin, AxesMax=amax, Color=c);
            saveFigures(this.figdir, closeFigure=true, prefix='call_apoe4_GAM_beta_apoe4');
    
            if opts.pvalues
                % PValue
                fprintf('NMFRadar.call_apoe4:\n')
                [h, crit_p, adj_ci_cvrg, Padj] = fdr_bh(apoe4.PValue', 0.05, 'dep', 'yes');
                fprintf("crit_p->%g; adj_ci_cvrg->%g", crit_p, adj_ci_cvrg)
                %disp(this.AxesLabels(h))
                disp(ascol(Padj(this.sorted_bases)))
                Padj = asrow(Padj); 
                amin = min(Padj, [], 'all');
                amax = max(Padj, [], 'all');
                if amin < 1e-15 || amax < 1e-15
                    return
                end
                figure
                axes_scaling = repmat({'log'}, [1 this.N_PATTERNS]);
                s2 = plot(this, Padj, AxesMin=amin, AxesMax=amax, ...
                    legend={'CDR=0, amy-'}, ti='FDR p-value \beta_{ApoE4}', AxesScaling=axes_scaling);
                %saveFigures(this.figdir, closeFigure=true, prefix='call_apoe4_FDR_p-value_beta_apoe4');
            end
        end

        function [s1,s2] = call_sex(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
                opts.pvalues logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            PC = ld.patternsofneurodegeneration20240630forimport;

            male = PC(contains(PC.ParamCoefficients, "sexM"), :);

            % Estimate
            P = male.Estimate';
            CI = 1.96*male.StdErr';
            amin = min(P-CI, [], 'all');
            amax = max(P+CI, [], 'all');
            figure
            c = this.suvr2viridis(min(P, [], "all"), suvr_min=-0.254, suvr_max=0.0108);
            s1 = plot_with_CI(this, P, CI, AxesMin=amin, AxesMax=amax, Color=c);
            saveFigures(this.figdir, closeFigure=true, prefix='call_sex_GAM_beta_sex');

            if opts.pvalues
                % PValue
                fprintf('NMFRadar.call_sex:\n')
                [h, crit_p, adj_ci_cvrg, Padj] = fdr_bh(male.PValue', 0.05, 'dep', 'yes');
                fprintf("crit_p->%g; adj_ci_cvrg->%g", crit_p, adj_ci_cvrg)
                %disp(this.AxesLabels(h))
                disp(ascol(Padj(this.sorted_bases)))
                Padj = asrow(Padj); 
                amin = min(Padj, [], 'all');
                amax = max(Padj, [], 'all');
                figure
                axes_scaling = repmat({'log'}, [1 this.N_PATTERNS]);
                s2 = plot(this, Padj, AxesMin=amin, AxesMax=amax, ...
                    legend={'CDR=0, amy-'}, ti='FDR p-value \beta_{sex}', AxesScaling=axes_scaling);
                %saveFigures(this.figdir, closeFigure=true, prefix='call_sex_FDR_p-value_beta_sex');
            end
        end    

        function [s1,s2] = call_groups(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
                opts.pvalues logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            CC = ld.patternsofneurodegeneration20240630forimport;

            assert(length(this.groups) == length(this.groupLabels))
            P = zeros(this.N_groups, this.N_PATTERNS);
            CI = zeros(this.N_groups, this.N_PATTERNS);
            h = nan(this.N_groups, this.N_PATTERNS);
            crit_p = nan(this.N_groups, 1);
            adj_ci_cvrg = nan(this.N_groups, 1);
            Padj = nan(this.N_groups, this.N_PATTERNS);
            selected = 2:length(this.groups);
            for ig = selected % cn is reference cohort/category
                try
                    T = CC(contains(CC.ParamCoefficients, this.groups{ig}), :);       
                    if opts.show
                        disp(T)
                    end

                    % Estimate
                    P(ig,:) = T.Estimate';
                    CI(ig,:) = 1.96*T.StdErr';
        
                    if opts.pvalues
                        % PValue
                        fprintf('NMFRadar.call_groups:\n')
                        [h__, crit_p(ig), adj_ci_cvrg(ig), Padj__] = ...
                            fdr_bh(T.PValue, 0.05, 'dep', 'yes');
                        fprintf("crit_p->%g; adj_ci_cvrg->%g", crit_p(ig), adj_ci_cvrg(ig))
                        h(ig,:) = asrow(h__);
                        Padj(ig,:) = asrow(Padj__);
                    end

                catch ME
                    handwarning(ME)
                end
            end
        
            h = figure;
            amin = min(P-CI, [], "all");
            amax = max(P+CI, [], "all");
            fprintf("%s: P: amin->%g, amax->%g\n", stackstr(), amin, amax)
            c = this.suvr2viridis(min(P(selected,:), [], 2), suvr_min=-0.254, suvr_max=0.0108);
            s1 = plot_groups_with_CI(this, P(selected,:), CI(selected,:), ...
                AxesMin=-0.3, AxesMax=0.05, ...
                AxesInterval=7, ...
                LineWidth=this.GROUP_LINEWIDTHS(selected), ...
                MinorGrid="off", MinorGridInterval=[], ...
                legend=this.groupLabels(selected), ...
                ti="", ...
                Color=c);
            saveFigure2(h, fullfile(this.figdir, 'call_groups_GAM_beta_groups'));

            if opts.pvalues
                % PValue
                U = table( ...
                    ascol(this.sorted_bases), ...
                    ascol(Padj(1,:)), ...
                    ascol(Padj(2,:)), ...
                    ascol(Padj(3,:)), ...
                    ascol(Padj(4,:)), ...
                    VariableNames={'sorted_bases', 'pval_preclin', 'pval_mci', 'pval_ad', 'pval_other'});
                disp(U)
                amin = min(Padj, [], "all");
                amax = max(Padj, [], "all");
                fprintf("%s: Padj: amin->%g, amax->%g\n", stackstr(), amin, amax)
                figure
                axes_scaling = repmat({'log'}, [1, this.N_PATTERNS]);
                s2 = plot_groups(this, Padj, ...
                    AxesMin=amin, AxesMax=amax, ...
                    legend=this.groupLabels(selected), ...
                    ti="FDR p-value \beta_groups", ...
                    AxesScaling=axes_scaling);
                saveFigures(this.figdir, closeFigure=true, prefix='call_groups_FDR_p-value_beta_groups');
            end
        end    

        function [s1,s2] = call_groups_percent(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
                opts.pvalues logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            CC = ld.patternsofneurodegeneration20240630forimport;

            % intercept
            intercept = CC(contains(CC.ParamCoefficients, "Intercept"), :);

            % groups
            assert(length(this.groups) == length(this.groupLabels))
            P = zeros(this.N_groups, this.N_PATTERNS);
            CI = zeros(this.N_groups, this.N_PATTERNS);
            h = nan(this.N_groups, this.N_PATTERNS);
            crit_p = nan(this.N_groups, 1);
            adj_ci_cvrg = nan(this.N_groups, 1);
            Padj = nan(this.N_groups, this.N_PATTERNS);
            selected = 2:length(this.groups);
            for ig = selected % cn is reference cohort/category
                try
                    T = CC(contains(CC.ParamCoefficients, this.groups{ig}), :);       
                    if opts.show
                        disp(T)
                    end

                    % Estimate percents with approx. of std error
                    R_ = T.Estimate'./(intercept.Estimate');
                    P(ig,:) = R_;
                    CI_ = R_.*sqrt((T.StdErr./T.Estimate).^2 + (intercept.StdErr./intercept.Estimate).^2)';
                    CI(ig,:) = 1.96*CI_';
        
                    if opts.pvalues
                        % PValue
                        fprintf('NMFRadar.call_groups:\n')
                        [h__, crit_p(ig), adj_ci_cvrg(ig), Padj__] = ...
                            fdr_bh(T.PValue, 0.05, 'dep', 'yes');
                        fprintf("crit_p->%g; adj_ci_cvrg->%g", crit_p(ig), adj_ci_cvrg(ig))
                        h(ig,:) = asrow(h__);
                        Padj(ig,:) = asrow(Padj__);
                    end

                catch ME
                    handwarning(ME)
                end
            end
        
            h = figure;
            amin = 100*min(P-CI, [], "all");
            amax = 100*max(P+CI, [], "all");
            fprintf("%s: P: amin->%g, amax->%g\n", stackstr(), amin, amax)
            c = this.suvr2viridis(min(100*P, [], 2), suvr_min=min(100*P, [], "all"), suvr_max=max(100*P, [], "all"));
            s1 = plot_groups_with_CI(this, 100*P(selected,:), 100*CI(selected,:), ...
                AxesMin=-25, AxesMax=5, ...
                AxesInterval=6, ...
                AxesPrecision=0, ...
                LineWidth=this.GROUP_LINEWIDTHS(selected), ...
                MinorGrid="off", MinorGridInterval=[], ...
                legend=this.groupLabels(selected), ...
                ti="", ...
                Color=c);
            saveFigure2(h, fullfile(this.figdir, 'call_groups_GAM_beta_groups_percent'));

            if opts.pvalues
                % PValue
                U = table( ...
                    ascol(this.sorted_bases), ...
                    ascol(Padj(1,:)), ...
                    ascol(Padj(2,:)), ...
                    ascol(Padj(3,:)), ...
                    ascol(Padj(4,:)), ...
                    VariableNames={'sorted_bases', 'pval_preclin', 'pval_mci', 'pval_ad', 'pval_other'});
                disp(U)
                amin = min(Padj, [], "all");
                amax = max(Padj, [], "all");
                fprintf("%s: Padj: amin->%g, amax->%g\n", stackstr(), amin, amax)
                figure
                axes_scaling = repmat({'log'}, [1, this.N_PATTERNS]);
                s2 = plot_groups(this, Padj, ...
                    AxesMin=amin, AxesMax=amax, ...
                    legend=this.groupLabels(selected), ...
                    ti="FDR p-value \beta_groups", ...
                    AxesScaling=axes_scaling);
                saveFigures(this.figdir, closeFigure=true, prefix='call_groups_FDR_p-value_beta_groups_percent');
            end
        end  

        function this = call_patt_weighted_fdg(this)
            c = 1;
            ld = load(fullfile(this.workdir, ...
                sprintf('baseline_%s', this.groups0{c}), ...
                sprintf('NumBases%i', this.N_bases_target), ...
                'components', this.matfile0));

            mu = nan(1, this.N_bases_target);
            sigma = nan(1, this.N_bases_target);
            for idx = 1:this.N_bases_target
                comp = ld.t.(sprintf("Components_%i", idx));
                mu(idx) = mean(comp, 1);
                sigma(idx) = std(comp, 1);
            end
            snr = mu./sigma;
            cov = sigma./mu;

            figure
            plot_with_CI(this, mu, sigma, ...
                AxesMin=dipmin(mu-sigma), AxesMax=dipmax(mu+sigma), ...
                legend=this.mergeDx(1), ...
                ti='Pattern-weighted FDG)');
            %saveFigures(this.figdir, closeFigure=true, prefix='patt_weighted_fdg_mu_sigma');

            figure
            plot(this, snr, ...
                AxesMin = dipmin(snr), AxesMax = dipmax(snr), ...
                legend=this.mergeDx(1), ...
                ti='\mu/\sigma Patterns of Covariance of FDG')
            %saveFigures(this.figdir, closeFigure=true, prefix='patt_weighted_fdg_snr');

            figure
            plot(this, cov, ...    
                AxesMin = dipmin(cov), AxesMax = dipmax(cov), ...
                legend=this.mergeDx(1), ...
                ti='\sigma/\mu Pattern of Covariance of FDG')
            %saveFigures(this.figdir, closeFigure=true, prefix='patt_weighted_fdg_cov');
        end
        
        function s = plot(this, P, opts)
            %% spider_plot_class_examples.m Example 5: Excel-like radar charts.

            arguments
                this mladni.NMFRadar %#ok<INUSA> 
                P double
                opts.legend {mustBeText} = {''}
                opts.ti {mustBeTextScalar} = ''
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesPrecision double = 3
                opts.AxesScaling cell = {}
                opts.AxesLabels = this.AxesLabels
                opts.Color double = []
                opts.LineWidth double = 2
            end
            if isempty(opts.Color)
                opts.Color = this.suvr2viridis(P);
            end
            P = P(this.sorted_bases);
            Nbases = size(P,2);
            
            % Delete variable in workspace if exists
            if exist('s', 'var')
                delete(s);
            end
            
            % Spider plot
            s = spider_plot_class(P);
            
            % Spider plot properties
            s.AxesLabels = opts.AxesLabels;
            s.AxesPrecision = opts.AxesPrecision;
            s.AxesDisplay = 'one';
            if ~isempty(opts.AxesMin) && ~isempty(opts.AxesMax)
                s.AxesLimits = [opts.AxesMin*ones(1,Nbases); opts.AxesMax*ones(1,Nbases)];
            end
            if ~isempty(opts.AxesInterval)
                s.AxesInterval = opts.AxesInterval;
            end
            %s.FillOption = 'on';
            %s.FillTransparency = 0.1;
            s.Color = opts.Color;
            s.LineWidth = opts.LineWidth;
            s.Marker = 'none';
            s.AxesFontSize = 11;
            s.LabelFontSize = 12;
            s.AxesColor = this.AXES_COLOR;
            s.AxesLabelsEdge = 'none';
            s.AxesLabelsRotate = 'on';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            if ~isemptytext(opts.legend)
                s.LegendLabels = opts.legend;
                s.LegendHandle.Location = 'northeastoutside';
                s.LegendHandle.FontSize = 13;
            end
            if ~isemptytext(opts.ti)
                title(opts.ti, 'FontSize', 14);
            end
        end

        function s = plot_with_CI(this, P, SE, opts)
            %% spider_plot_class_examples.m 
            %  Example 5 with Excel-like radar charts and
            %  Example 9 with shaded areas

            arguments
                this mladni.NMFRadar %#ok<INUSA> 
                P double
                SE double
                opts.legend {mustBeText} = {''}
                opts.ti {mustBeTextScalar} = ''
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesPrecision double = 2
                opts.AxesScaling cell = {}
                opts.AxesLabels = this.AxesLabelsNull
                opts.Color double = []
                opts.LineWidth double = 1
            end
            if isempty(opts.Color)
                opts.Color = this.suvr2viridis(P);
            end
            P = P(this.sorted_bases);
            SE = SE(this.sorted_bases);

            axes_shaded_limits = { ...
                [asrow(P) - asrow(SE); asrow(P) + asrow(SE)]};
            Nbases = size(P,2);
            
            % Delete variable in workspace if exists
            if exist('s', 'var')
                delete(s);
            end
            
            % Spider plot
            try
                s = spider_plot_class(P);
            catch ME
                handwarning(ME)
            end
            
            % Spider plot properties
            if ~isempty(opts.AxesMin) && ~isempty(opts.AxesMax)
                s.AxesLimits = [opts.AxesMin*ones(1,Nbases); opts.AxesMax*ones(1,Nbases)];
            end
            if ~isempty(opts.AxesInterval)
                s.AxesInterval = opts.AxesInterval;
            end            
            s.AxesLabels = opts.AxesLabels;
            s.AxesShaded = 'on';
            s.AxesShadedLimits = axes_shaded_limits;
            s.AxesShadedColor = {opts.Color};
            s.AxesShadedTransparency = this.ERR_ALPHA;
            s.AxesPrecision = opts.AxesPrecision;
            s.AxesDisplay = 'one';
            %s.FillOption = 'on';
            %s.FillTransparency = 0.1;
            s.Color = opts.Color;
            s.LineWidth = opts.LineWidth;
            s.Marker = 'none';
            s.AxesFontSize = 11;
            s.LabelFontSize = 12;
            s.AxesColor = this.AXES_COLOR;
            s.AxesLabelsEdge = 'none';
            s.AxesLabelsRotate = 'on';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            if ~isemptytext(opts.legend)
                s.LegendLabels = opts.legend;
                s.LegendHandle.Location = 'northeastoutside';
                s.LegendHandle.FontSize = 13;
            end
            if ~isemptytext(opts.ti)
                title(opts.ti, 'FontSize', 14);
            end
        end

        function s = plot_groups(this, P, opts)
            %% spider_plot_class_examples.m Example 5: Excel-like radar charts.

            arguments
                this mladni.NMFRadar %#ok<INUSA> 
                P double
                opts.legend {mustBeText} = {''}
                opts.ti {mustBeTextScalar} = ''
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesPrecision double = 3
                opts.AxesScaling cell = {}
                opts.Color double = []
                opts.LineStyle {mustBeText} = '-'
                opts.LineWidth double = 1
                opts.AxesLabels = this.AxesLabels
                opts.MinorGrid = "on"
                opts.MinorGridInterval double = [];
            end
            if isempty(opts.Color)
                opts.Color = this.suvr2viridis(P);
            end
            P = P(:, this.sorted_bases);
            Nbases = size(P,2);
            
            % Delete variable in workspace if exists
            if exist('s', 'var')
                delete(s);
            end
            
            % Spider plot
            s = spider_plot_class(P);
            
            % Spider plot properties
            if ~isempty(opts.AxesMin) && ~isempty(opts.AxesMax)
                s.AxesLimits = [opts.AxesMin*ones(1,Nbases); opts.AxesMax*ones(1,Nbases)];
            end
            if ~isempty(opts.AxesInterval)
                s.AxesInterval = opts.AxesInterval;
            end
            s.AxesLabels = opts.AxesLabels;
            s.AxesPrecision = opts.AxesPrecision;
            s.AxesDisplay = 'one';
            %s.FillOption = 'on';
            %s.FillTransparency = 0.1;
            s.Color = opts.Color;
            s.LineStyle = opts.LineStyle;
            s.LineWidth = opts.LineWidth;
            s.Marker = 'none';
            s.AxesFontSize = 11;
            s.LabelFontSize = 12;
            s.AxesColor = this.AXES_COLOR;
            s.AxesLabelsEdge = 'none';
            s.AxesLabelsRotate = 'on';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            if ~isemptytext(opts.legend)
                s.LegendLabels = opts.legend;
                s.LegendHandle.Location = 'northeastoutside';
                s.LegendHandle.FontSize = 13;
            end
            if ~isemptytext(opts.ti)
                title(opts.ti, 'FontSize', 14);
            end
        end
        function s = plot_groups_with_CI(this, P, SE, opts)
            %% spider_plot_class_examples.m 
            %  Example 5 with Excel-like radar charts and
            %  Example 9 with shaded areas

            arguments
                this mladni.NMFRadar %#ok<INUSA> 
                P double
                SE double
                opts.legend {mustBeText} = {''}
                opts.ti {mustBeTextScalar} = ''
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesPrecision double = 2
                opts.AxesScaling cell = {}
                opts.Color double = []
                opts.LineStyle {mustBeText} = '-'
                opts.LineWidth double = 1
                opts.AxesLabels = this.AxesLabelsNull
                opts.MinorGrid = "on"
                opts.MinorGridInterval double = [];
            end
            if isempty(opts.Color)
                opts.Color = this.suvr2viridis(P);
            end
            P = P(:, this.sorted_bases);
            SE = SE(:, this.sorted_bases);
            axes_shaded_limits = { ...
                [P(1,:) - SE(1,:); P(1,:) + SE(1,:)], ...
                [P(2,:) - SE(2,:); P(2,:) + SE(2,:)], ...
                [P(3,:) - SE(3,:); P(3,:) + SE(3,:)]};  % KLUDGE
            axes_shaded_colors = { ...
                opts.Color(1,:), opts.Color(2,:), opts.Color(3,:)};  % KLUDGE
            Nbases = size(P,2);
            
            % Delete variable in workspace if exists
            if exist('s', 'var')
                delete(s);
            end
            
            % Spider plot
            try
                s = spider_plot_class(P);
            catch ME
                handwarning(ME)
            end
            
            % Spider plot properties
            if ~isempty(opts.AxesMin) && ~isempty(opts.AxesMax)
                s.AxesLimits = [opts.AxesMin*ones(1,Nbases); opts.AxesMax*ones(1,Nbases)];
            end
            if ~isempty(opts.AxesInterval)
                s.AxesInterval = opts.AxesInterval;
            end
            s.MinorGrid = opts.MinorGrid;
            if ~isempty(opts.MinorGridInterval)
                s.MinorGridInterval = opts.MinorGridInterval;
            end
            s.AxesLabels = opts.AxesLabels;
            s.AxesShaded = 'on';
            s.AxesShadedLimits = axes_shaded_limits;
            s.AxesShadedColor = axes_shaded_colors;
            s.AxesShadedTransparency = this.ERR_ALPHA;
            s.AxesPrecision = opts.AxesPrecision;
            s.AxesDisplay = 'one';
            %s.FillOption = 'on';
            %s.FillTransparency = 0.1;
            s.Color = opts.Color;
            s.LineStyle = opts.LineStyle;
            s.LineWidth = opts.LineWidth;
            s.Marker = 'none';
            s.AxesFontSize = 11;
            s.LabelFontSize = 12;
            s.AxesColor = this.AXES_COLOR;
            s.AxesLabelsEdge = 'none';
            s.AxesLabelsRotate = 'on';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            if ~isemptytext(opts.legend)
                s.LegendLabels = opts.legend;
                s.LegendHandle.Location = 'northeastoutside';
                s.LegendHandle.FontSize = 13;
            end
            if ~isemptytext(opts.ti)
                title(opts.ti, 'FontSize', 14);
            end
        end
        function h = plot_beta0_to_beta1(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show_labels logical = true
            end

            pwd0 = pushd(this.workdir);

            %% unwraps the radar plots to show panels of maps of beta0 -> beta1_{Dx group}
            NP = this.N_PATTERNS;
            ld = load(this.matfile_cohort);
            CC = ld.patternsofneurodegeneration20240630forimport;
            CC = natsortrows(CC);
            CC(contains(CC.ParamCoefficients, "apoe4"),:) = [];
            CC(contains(CC.ParamCoefficients, "sexM"),:) = [];
            % CC(contains(CC.ParamCoefficients, "cohortCDR>0,amy-"),:) = [];
            indices = nan(1, NP);
            %for idx = 1:NP
            %    indices(this.sorted_bases(idx)) = idx;
            %end
            indices(this.sorted_bases) = 1:NP;
            indices = num2cell(indices);
            if opts.show_labels
                labels = cellfun(@(x) sprintf('p%i', x), indices, UniformOutput=false);
                filesuffix = '_labelled';
            else
                labels = cellfun(@(x) '', indices);
                filesuffix = '';
            end            
            position = [1,1,1045,1692];
            cmin = -0.254;
            cmax = 0.0108;
            marker_size = 100;
            fscale = 3;

            figure        
            scatter(CC.Estimate(1:NP), CC.Estimate(  NP+1:2*NP), marker_size, CC.Estimate(  NP+1:2*NP), "filled");
            colormap(viridis);
            clim([cmin, cmax]);
            % labelpoints(CC.Estimate(1:NP), CC.Estimate(  NP+1:2*NP), labels, 'E')
            xlim([0.7 1.4])
            ylim([-0.3 0.02])
            % xlabel("\beta_{CDR=0,amy-}", FontSize=14)
            % ylabel("\beta_{CDR=0,amy+}", FontSize=14)
            fontsize(gcf, scale=fscale)
            grid on
            pbaspect([1 3 1])
            set(gcf, position=position)
            saveFigure2(gcf, fullfile(this.figdir, "plot_beta0_to_beta1_cdr=0_amy+" + filesuffix))

            figure        
            scatter(CC.Estimate(1:NP), CC.Estimate(3*NP+1:4*NP), marker_size, CC.Estimate(3*NP+1:4*NP), "filled")    
            colormap(viridis);   
            clim([cmin, cmax]);     
            % labelpoints(CC.Estimate(1:NP), CC.Estimate(3*NP+1:4*NP), labels, 'E')
            xlim([0.7 1.4])
            ylim([-0.3 0.02])
            % yticklabels([])
            % xlabel("\beta_{CDR=0,amy-}", FontSize=14)
            % ylabel("\beta_{CDR>0,amy-}", FontSize=14)
            fontsize(gcf, scale=fscale)
            grid on
            pbaspect([1 3 1])
            set(gcf, position=position)
            saveFigure2(gcf, fullfile(this.figdir, "plot_beta0_to_beta1_cdr_gt_0_amy-" + filesuffix))

            figure
            scatter(CC.Estimate(1:NP), CC.Estimate(2*NP+1:3*NP), marker_size, CC.Estimate(2*NP+1:3*NP), "filled")
            colormap(viridis);
            clim([cmin, cmax]);
            % labelpoints(CC.Estimate(1:NP), CC.Estimate(2*NP+1:3*NP), labels, 'E')
            xlim([0.7 1.4])
            ylim([-0.3 0.02])
            % yticklabels([])
            % xlabel("\beta_{CDR=0,amy-}", FontSize=14)
            % ylabel("\beta_{CDR=0.5,amy+}", FontSize=14)
            fontsize(gcf, scale=fscale)
            grid on
            pbaspect([1 3 1])
            set(gcf, position=position)
            saveFigure2(gcf, fullfile(this.figdir, "plot_beta0_to_beta1_cdr=0p5_amy+" + filesuffix))

            figure
            scatter(CC.Estimate(1:NP), CC.Estimate(4*NP+1:5*NP), marker_size, CC.Estimate(4*NP+1:5*NP), "filled")
            colormap(viridis);
            clim([cmin, cmax]);
            % labelpoints(CC.Estimate(1:NP), CC.Estimate(4*NP+1:5*NP), labels, 'E')
            xlim([0.7 1.4])
            ylim([-0.3 0.02])
            % yticklabels([])
            % xlabel("\beta_{CDR=0,amy-}", FontSize=14)
            % ylabel("\beta_{CDR>0.5,amy+}", FontSize=14)
            fontsize(gcf, scale=fscale)
            grid on
            pbaspect([1 3 1])
            set(gcf, position=position)
            saveFigure2(gcf, fullfile(this.figdir, "plot_beta0_to_beta1_cdr_gt_0p5_amy+" + filesuffix))

            popd(pwd0);
        end

        function c = suvr2viridis(this, suvr, opts)
            arguments
                this mladni.NMFRadar
                suvr double = 0
                opts.suvr_max = []
                opts.suvr_min = []
            end
            if isempty(opts.suvr_max)
                opts.suvr_max = max(suvr, [], "all");
            end
            if isempty(opts.suvr_min)
                opts.suvr_min = min(suvr, [], "all");
            end

            v = viridis(1000);
            suvr_range = opts.suvr_max - opts.suvr_min;
            suvr_frac = (suvr - opts.suvr_min)/suvr_range;
            idx = round((1 - suvr_frac) + 1000*suvr_frac);
            idx(idx < 1) = 1;
            idx(idx > 1000) = 1000;
            c = v(idx, :);
        end

        function T = table_cohort_coefficients_sorted(this)
            %% sorts this.matfile_cohort and writes this.matfile_cohort_sorted, along with corresponding csv.

            patterns_are_adjacent = false;

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            T0 = ld.patternsofneurodegeneration20240630forimport;
            T = T0;

            N_bundle = size(T0, 1) / 24;  % bundle of covariates:  Intercept, SexM, apoe4, chohortCDR*, ...
            assert(rem(size(T0, 1), 24) == 0)
            if patterns_are_adjacent
                for b = 1:N_bundle
                    rows = (b - 1)*24 + (1:24);
                    U0 = T0(rows, :);
                    U1 = U0(this.sorted_bases, :);
                    U1.ParamCoefficients = U0.ParamCoefficients;
                    T(rows,:) = U1;
                end
            else
                for b = 1:N_bundle
                    rows = ascol(b + (0:23)*N_bundle);
                    U0 = T0(rows, :);
                    U1 = U0(this.sorted_bases, :);
                    U1.ParamCoefficients = U0.ParamCoefficients;  % sort the Estimates, then replace the ParamCoeff labels
                    T(rows,:) = U1;
                end
            end

            % FDR by Benjamini-Hochberg
            [~,~,~,FdrPValue] = fdr_bh(T.PValue, 0.05, 'dep', 'yes');
            T = addvars(T, FdrPValue);
            T = natsortrows(T);
            
            save(fullfile(this.workdir, this.matfile_cohort_sorted), "T");
            writetable(T, fullfile(this.workdir, myfileprefix(this.matfile_cohort_sorted) + ".csv"));
        end
        function T = table_patt_weighted_fdg(this)

            nmfh = mladni.NMFHierarchies();
            T = nmfh.table_patt_weighted_fdg();
        end

        function this = NMFRadar(varargin)
            this.workdir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
            
            if ~isfile(fullfile(this.workdir, this.matfile_cohort_sorted))
                this.table_cohort_coefficients_sorted();
            end
        end
    end

    %% PROTECTED
    
    properties (Access = protected)
        sorted_bases_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
