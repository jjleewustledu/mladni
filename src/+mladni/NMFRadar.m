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
        N_PATTERNS = mladni.NMF.N_PATTERNS
    end

    properties
        groups = { ...
            'CDR=0,amy+' 'CDR=0.5,amy+' 'CDR>0.5,amy+'} 
        groupLabels = {...
            'CDR=0,amy+' 'CDR=0.5,amy+' 'CDR>0.5,amy+'}
        groups0 = { ...
            'cn', 'preclinical', ...
            'cdr_0p5_apos', ...
            'cdr_gt_0p5_apos', ...
            'cdr_gt_0_aneg'}
        mergeDx = { ...
            'CDR=0,amy-' 'CDR=0,amy+' 'CDR=0.5,amy+' 'CDR>0.5,amy+'}
        
        matfile0 = 'NMFCovariates_table_covariates_1stscan_longitudinal.mat'
        matfile_cohort = 'CohortCoefficients_20230928.mat' 
        % Output from R: Singularity/ADNI/NMF_FDG/patterns_of_neurodegeneration_20230921.Rmd
        % b2 <- gam(list(
        % y0~s(age,k=20,by=interaction(sex))+sex+apoe4+cohort, y1~s(age,k=20,by=interaction(sex))+sex+apoe4+cohort, y2~s(age,k=20,by=interaction(sex))+sex+apoe4+cohort, ...
        % family=mvn(d=24), data=soto)
        % b2
        % summary(b2)  

        workdir
    end

    properties (Dependent)
        AxesLabels
        AxesLabelsNull
        figdir
        label_permute
        N_bases_target
        N_groups
        sorted_bases % ordered by pattern-weighted FDG SUvR
    end

    methods % GET
        function g = get.AxesLabelsNull(~)
            g = {'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''};
        end
        function g = get.AxesLabels(this)
            g = {'P1', 'P2', 'P3', 'P4', ...
                 'P5', 'P6', 'P7', 'P8', ...
                 'P9', 'P10', 'P11', 'P12', ...
                 'P13', 'P14', 'P15', 'P16', ...
                 'P17', 'P18', 'P19', 'P20', ...
                 'P21', 'P22', 'P23', 'P24'};
            g = g(this.sorted_bases);
        end
        function g = get.sorted_bases(this)
            if ~isempty(this.sorted_bases_)
                g = this.sorted_bases_;
                return
            end
            T = this.table_patt_weighted_fdg();
            g = asrow(T.sorted_bases);
        end
        function g = get.figdir(this)
            g = fullfile(this.workdir, 'baseline_cn', 'results');
        end
        function g = get.label_permute(this)
            switch this.N_PATTERNS
                case 16
                    g = [1 2 13 16 3 4 5 6 7 8 9 10 11 12 14 15]; % ParamCoeff indices are 0, 1, 10, 11, 12, ...                    
                case 22
                    g = [1 2 13 16 17 18 19 20 21 22 3 4 5 6 7 8 9 10 11 12 14 15]; % ParamCoeff indices are 0, 1, 10, 11, 12, ...
                case 24
                    g = [1 2 13 16 17 18 19 20 21 22 23 24 3 4 5 6 7 8 9 10 11 12 14 15]; % ParamCoeff indices are 0, 1, 10, 11, 12, ...
                otherwise
                    error("mladni:ValueError", "%s: N_PATTERNS->%i", stackstr(), mladni.NMF.N_PATTERNS);
            end
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
                opts.pvalues logical = true
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            PC = ld.CohortCoefficients20230928;
            intercept = PC(contains(PC.ParamCoefficients, "Intercept"), :);
            if opts.show
                disp(intercept)
            end

            % Estimate
            P = intercept.Estimate';
            SE = intercept.StdErr';
            amin = min(P-SE, [], 'all');
            amax = max(P+SE, [], 'all');
            figure
            s1 = plot_with_stderr(this, P, SE, AxesMin=amin, AxesMax=amax, ...
                legend={'CDR=0, amy-'}, ti='');
            saveFigures(this.figdir, closeFigure=true, prefix='GAM beta_intercept');
    
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
                saveFigures(this.figdir, closeFigure=true, prefix='FDR p-value beta_intercept');
            end
        end
        function [s1,s2] = call_apoe4(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
                opts.pvalues logical = true
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            PC = ld.CohortCoefficients20230928;
            apoe4 = PC(contains(PC.ParamCoefficients, "apoe4"), :);
            if opts.show
                disp(apoe4)
            end

            % Estimate
            %load("P.mat");
            %load("SE.mat");
            P = apoe4.Estimate';
            SE = apoe4.StdErr';
            amin = min(P-SE, [], 'all');
            amax = max(P+SE, [], 'all');
            figure
            s1 = plot_with_stderr(this, P, SE, AxesMin=amin, AxesMax=amax, ...
                legend={'CDR=0, amy-'}, ti='');
            saveFigures(this.figdir, closeFigure=true, prefix='GAM beta_apoe4');
    
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
                saveFigures(this.figdir, closeFigure=true, prefix='FDR p-value beta_apoe4');
            end
        end
        function [s1,s2] = call_sex(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
                opts.pvalues logical = true
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            PC = ld.CohortCoefficients20230928;

            male = PC(contains(PC.ParamCoefficients, "sexM"), :);
            %male = sortrows(male, "ParamCoefficients");
            %male = male(this.label_permute, :);

            if opts.show
                disp(male)
            end

            % Estimate
            P = male.Estimate';
            SE = male.StdErr';
            amin = min(P-SE, [], 'all');
            amax = max(P+SE, [], 'all');
            figure
            s1 = plot_with_stderr(this, P, SE, AxesMin=amin, AxesMax=amax, ...
                legend={'CDR=0, amy-'}, ti='');
            saveFigures(this.figdir, closeFigure=true, prefix='GAM beta_sex');

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
                saveFigures(this.figdir, closeFigure=true, prefix='FDR p-value beta_sex');
            end
        end      
        function [s1,s2] = call_groups(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
                opts.pvalues logical = true
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            CC = ld.CohortCoefficients20230928;

            assert(length(this.groups) == length(this.groupLabels))
            P = zeros(this.N_groups+1, this.N_PATTERNS);
            SE = zeros(this.N_groups+1, this.N_PATTERNS);
            h = nan(this.N_groups, this.N_PATTERNS);
            crit_p = nan(this.N_groups, 1);
            adj_ci_cvrg = nan(this.N_groups, 1);
            Padj = nan(this.N_groups, this.N_PATTERNS);
            for ig = 1:length(this.groups) % cn is reference cohort/category
                try
                    T = CC(contains(CC.ParamCoefficients, this.groups{ig}), :);
                    %T = sortrows(T, "ParamCoefficients");
                    %T = T(this.label_permute, :);        
                    if opts.show
                        disp(T)
                    end

                    % Estimate
                    P(ig,:) = T.Estimate';
                    SE(ig,:) = T.StdErr';
        
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
        
            figure
            amin = min(P-SE, [], "all");
            amax = max(P+SE, [], "all");
            fprintf("%s: P: amin->%g, amax->%g\n", stackstr(), amin, amax)
            s1 = plot_groups_with_stderr(this, P([1,4,2,3],:), SE([1,4,2,3],:), ...
                AxesMin=-0.3, AxesMax=0.05, ...
                AxesInterval=7, ...
                Color=[0.85,0.325,0.098;0,0,0;0.929,0.674,0.125;0.494,0.184,0.556], ...
                LineWidth=[2,0.5,2,2], ...
                MinorGrid="off", MinorGridInterval=[], ...
                legend=[this.groupLabels{1}, '\beta=0',this.groupLabels(2:3)], ...
                ti="");
            saveFigures(this.figdir, closeFigure=true, prefix='GAM beta_groups');

            if opts.pvalues
                % PValue
                U = table( ...
                    ascol(this.sorted_bases), ...
                    ascol(Padj(1,:)), ...
                    ascol(Padj(2,:)), ...
                    ascol(Padj(3,:)), ...
                    VariableNames={'sorted_bases', 'pval_preclin', 'pval_mci', 'pval_ad'});
                disp(U)
                amin = min(Padj, [], "all");
                amax = max(Padj, [], "all");
                fprintf("%s: Padj: amin->%g, amax->%g\n", stackstr(), amin, amax)
                figure
                axes_scaling = repmat({'log'}, [1, this.N_PATTERNS]);
                s2 = plot_groups(this, Padj, ...
                    AxesMin=amin, AxesMax=amax, ...
                    legend=this.groupLabels, ...
                    ti="FDR p-value \beta_groups", ...
                    AxesScaling=axes_scaling);
                saveFigures(this.figdir, closeFigure=true, prefix='FDR p-value beta_groups');
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
            plot_with_stderr(this, mu, sigma, ...
                AxesMin=dipmin(mu-sigma), AxesMax=dipmax(mu+sigma), ...
                legend=this.mergeDx(1), ...
                ti='Pattern-weighted FDG (SUVR)');
            saveFigures(this.figdir, closeFigure=true, prefix='patt_weighted_fdg_mu_sigma');

            figure
            plot(this, snr, ...
                AxesMin = dipmin(snr), AxesMax = dipmax(snr), ...
                legend=this.mergeDx(1), ...
                ti='\mu/\sigma Pattern-Weighted FDG (SUVR)')
            saveFigures(this.figdir, closeFigure=true, prefix='patt_weighted_fdg_snr');

            figure
            plot(this, cov, ...    
                AxesMin = dipmin(cov), AxesMax = dipmax(cov), ...
                legend=this.mergeDx(1), ...
                ti='\sigma/\mu Pattern-Weighted FDG (SUVR)')
            saveFigures(this.figdir, closeFigure=true, prefix='patt_weighted_fdg_cov');
        end
        
        function s = plot(this, P, opts)
            %% spider_plot_class_examples.m Example 5: Excel-like radar charts.

            arguments
                this mladni.NMFRadar %#ok<INUSA> 
                P double
                opts.legend {mustBeText} = this.mergeDx(1)
                opts.ti {mustBeTextScalar}
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesScaling cell = {}
                opts.AxesLabels = this.AxesLabels
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
            s.AxesPrecision = 3;
            s.AxesDisplay = 'one';
            if ~isempty(opts.AxesMin) && ~isempty(opts.AxesMax)
                s.AxesLimits = [opts.AxesMin*ones(1,Nbases); opts.AxesMax*ones(1,Nbases)];
            end
            if ~isempty(opts.AxesInterval)
                s.AxesInterval = opts.AxesInterval;
            end
            %s.FillOption = 'on';
            %s.FillTransparency = 0.1;
            s.Color = [139, 0, 0; 240, 128, 128]/255;
            s.LineWidth = 2;
            s.Marker = 'none';
            s.AxesFontSize = 11;
            s.LabelFontSize = 12;
            s.AxesColor = [0.8, 0.8, 0.8];
            s.AxesLabelsEdge = 'none';
            s.AxesLabelsRotate = 'on';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            s.LegendLabels = opts.legend;
            s.LegendHandle.Location = 'northeastoutside';
            s.LegendHandle.FontSize = 13;
            title(opts.ti, 'FontSize', 14);
        end
        function s = plot_with_stderr(this, P, SE, opts)
            %% spider_plot_class_examples.m 
            %  Example 5 with Excel-like radar charts and
            %  Example 9 with shaded areas

            arguments
                this mladni.NMFRadar %#ok<INUSA> 
                P double
                SE double
                opts.legend {mustBeText} = this.mergeDx(1)
                opts.ti {mustBeTextScalar}
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesScaling cell = {}
                opts.AxesLabels = this.AxesLabelsNull
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
            s.AxesShadedColor = 'b';
            s.AxesShadedTransparency = 0.1;
            s.AxesPrecision = 2;
            s.AxesDisplay = 'one';
            %s.FillOption = 'on';
            %s.FillTransparency = 0.1;
            %s.Color = [0, 0, 139; 128, 128, 240]/255;
            s.LineWidth = 2;
            s.Marker = 'none';
            s.AxesFontSize = 11;
            s.LabelFontSize = 12;
            s.AxesColor = [0.8, 0.8, 0.8];
            s.AxesLabelsEdge = 'none';
            s.AxesLabelsRotate = 'on';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            s.LegendLabels = opts.legend;
            s.LegendHandle.Location = 'northeastoutside';
            s.LegendHandle.FontSize = 13;
            title(opts.ti, 'FontSize', 14);
        end
        function s = plot_groups(this, P, opts)
            %% spider_plot_class_examples.m Example 5: Excel-like radar charts.

            arguments
                this mladni.NMFRadar %#ok<INUSA> 
                P double
                opts.legend {mustBeText} = this.mergeDx(2:end)
                opts.ti {mustBeTextScalar}
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesScaling cell = {}
                opts.Color double = []
                opts.LineStyle {mustBeText} = '-'
                opts.LineWidth double = 2
                opts.AxesLabels = this.AxesLabels
                opts.MinorGrid = "on"
                opts.MinorGridInterval double = [];
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
            if ~isempty(opts.Color)
                s.Color = opts.Color;
            end
            if ~isempty(opts.LineStyle)
                s.LineStyle = opts.LineStyle;
            end
            if ~isempty(opts.AxesMin) && ~isempty(opts.AxesMax)
                s.AxesLimits = [opts.AxesMin*ones(1,Nbases); opts.AxesMax*ones(1,Nbases)];
            end
            if ~isempty(opts.AxesInterval)
                s.AxesInterval = opts.AxesInterval;
            end
            s.AxesLabels = opts.AxesLabels;
            s.AxesPrecision = 3;
            s.AxesDisplay = 'one';
            %s.FillOption = 'on';
            %s.FillTransparency = 0.1;
            s.Color = [139, 0, 0; 240, 128, 128]/255;
            s.LineWidth = 2;
            s.Marker = 'none';
            s.AxesFontSize = 11;
            s.LabelFontSize = 12;
            s.AxesColor = [0.8, 0.8, 0.8];
            s.AxesLabelsEdge = 'none';
            s.AxesLabelsRotate = 'on';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            s.LegendLabels = opts.legend;
            s.LegendHandle.Location = 'northeastoutside';
            s.LegendHandle.FontSize = 13;
            title(opts.ti, 'FontSize', 14);
        end
        function s = plot_groups_with_stderr(this, P, SE, opts)
            %% spider_plot_class_examples.m 
            %  Example 5 with Excel-like radar charts and
            %  Example 9 with shaded areas

            arguments
                this mladni.NMFRadar %#ok<INUSA> 
                P double
                SE double
                opts.legend {mustBeText} = this.mergeDx(2:end)
                opts.ti {mustBeTextScalar}
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesScaling cell = {}
                opts.Color double = []
                opts.LineStyle {mustBeText} = '-'
                opts.LineWidth double = 2
                opts.AxesLabels = this.AxesLabelsNull
                opts.MinorGrid = "on"
                opts.MinorGridInterval double = [];
            end
            P = P(:, this.sorted_bases);
            SE = SE(:, this.sorted_bases);

            axes_shaded_limits = { ...
                [P(1,:) - SE(1,:); P(1,:) + SE(1,:)], ...
                [P(2,:) - SE(2,:); P(2,:) + SE(2,:)], ...
                [P(3,:) - SE(3,:); P(3,:) + SE(3,:)]};
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
            if ~isempty(opts.Color)
                s.Color = opts.Color;
            end
            if ~isempty(opts.LineStyle)
                s.LineStyle = opts.LineStyle;
            end
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
            s.AxesShadedColor = 'b';
            s.AxesShadedTransparency = 0.1;
            s.AxesPrecision = 2;
            s.AxesDisplay = 'one';
            %s.FillOption = 'on';
            %s.FillTransparency = 0.1;
            %s.Color = [0, 0, 139; 128, 128, 240]/255;
            s.LineWidth = opts.LineWidth;
            s.Marker = 'none';
            s.AxesFontSize = 11;
            s.LabelFontSize = 12;
            s.AxesColor = [0.8, 0.8, 0.8];
            s.AxesLabelsEdge = 'none';
            s.AxesLabelsRotate = 'on';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            s.LegendLabels = opts.legend;
            s.LegendHandle.Location = 'northeastoutside';
            s.LegendHandle.FontSize = 13;
            title(opts.ti, 'FontSize', 14);
        end
        function h = plot_beta0_to_beta1(this)
            NP = this.N_PATTERNS;
            ld = load(this.matfile_cohort);
            CC = ld.CohortCoefficients20230928;
            CC(contains(CC.ParamCoefficients, "apoe4"),:) = [];
            CC(contains(CC.ParamCoefficients, "sexM"),:) = [];
            CC(contains(CC.ParamCoefficients, "cohortCDR>0,amy-"),:) = [];
            indices = nan(1, NP);
            %for idx = 1:NP
            %    indices(this.sorted_bases(idx)) = idx;
            %end
            indices(this.sorted_bases) = 1:NP;
            indices = num2cell(indices);
            labels = cellfun(@(x) sprintf('P%i', x), indices, UniformOutput=false);

            figure        
            plot(CC.Estimate(1:NP), CC.Estimate(  NP+1:2*NP), LineStyle="none", Marker="o", MarkerSize=24)            
            labelpoints(CC.Estimate(1:NP), CC.Estimate(  NP+1:2*NP), labels, 'C', 0.25, 1)
            xlim([0.7 1.4])
            ylim([-0.3 0.015])
            xlabel("\beta_{CDR=0,amy-}", FontSize=14)
            ylabel("\beta_{CDR=0,amy+}", FontSize=14)
            fontsize(gcf, scale=1.2)
            grid on
            pbaspect([1 1.618 1])
            set(gcf, position=[100,100,618,1000])
            saveFigure2(gcf, fullfile(this.figdir, "preclinical"))

            figure
            plot(CC.Estimate(1:NP), CC.Estimate(2*NP+1:3*NP), LineStyle="none", Marker="o", MarkerSize=24)
            labelpoints(CC.Estimate(1:NP), CC.Estimate(2*NP+1:3*NP), labels, 'C', 0.25, 1)
            xlim([0.7 1.4])
            ylim([-0.3 0.015])
            xlabel("\beta_{CDR=0,amy-}", FontSize=14)
            ylabel("\beta_{CDR=0.5,amy+}", FontSize=14)
            fontsize(gcf, scale=1.2)
            grid on
            pbaspect([1 1.618 1])
            set(gcf, position=[100,100,618,1000])
            saveFigure2(gcf, fullfile(this.figdir, "mci"))

            figure
            plot(CC.Estimate(1:NP), CC.Estimate(3*NP+1:4*NP), LineStyle="none", Marker="o", MarkerSize=24)
            labelpoints(CC.Estimate(1:NP), CC.Estimate(3*NP+1:4*NP), labels, 'C', 0.25, 1)
            xlim([0.7 1.4])
            ylim([-0.3 0.015])
            xlabel("\beta_{CDR=0,amy-}", FontSize=14)
            ylabel("\beta_{CDR>0.5,amy+}", FontSize=14)
            fontsize(gcf, scale=1.2)
            grid on
            pbaspect([1 1.618 1])
            set(gcf, position=[100,100,618,1000])
            saveFigure2(gcf, fullfile(this.figdir, "ad"))
        end
        
        function h = plot_beta0_for_groups(this)
        end
        function h = plot_beta1_for_groups(this)
        end

        function T = table_patt_weighted_fdg(this)
            c = 1;
            ld = load(fullfile(this.workdir, ...
                sprintf('baseline_%s', this.groups0{c}), ...
                sprintf('NumBases%i', this.N_bases_target), ...
                'components', this.matfile0));

            indices_bases = 1:this.N_bases_target; %#ok<PROP>
            mu = nan(1, this.N_bases_target);
            sigma = nan(1, this.N_bases_target);
            for idx = 1:this.N_bases_target
                comp = ld.t.(sprintf("Components_%i", idx));
                mu(idx) = mean(comp, 1);
                sigma(idx) = std(comp, 1);
            end

            T = table(indices_bases', mu', sigma', VariableNames={'indices_bases', 'mu', 'sigma'}); %#ok<PROP>
            T = sortrows(T, 2, "descend");
            this.sorted_bases_ = T.indices_bases;
        end

        function this = NMFRadar(varargin)
            this.workdir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
        end
    end

    %% PROTECTED
    
    properties (Access = protected)
        sorted_bases_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
