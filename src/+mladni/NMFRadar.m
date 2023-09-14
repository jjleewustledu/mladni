classdef NMFRadar
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
            'CDR=0,amy-' 'CDR=0,amy+' 'CDR=0.5,amy+' 'CDR>0.5,amy+' 'CDR>0,amy-'} % DEPRECATED
        groupLabels = {...
            'CDR=0,amy-' 'CDR=0,amy+' 'CDR=0.5,amy+' 'CDR>0.5,amy+' 'CDR>0,amy-'} % DEPRECATED
        groups0 = { ...
            'cn', 'preclinical', ...
            'cdr_0p5_apos', ...
            'cdr_gt_0p5_apos', ...
            'cdr_gt_0_aneg'}
        mergeDx = { ...
            'CDR=0,amy-' 'CDR=0,amy+' 'CDR=0.5,amy+' 'CDR>0.5,amy+' 'CDR>0,amy-'}
        
        matfile0 = 'NMFCovariates_table_covariates_1stscan_longitudinal.mat'
        matfile_cohort = 'CohortCoefficients_20230721.mat' % output from R:patterns_of_neurodegeneration*.Rmd
        matfile_param = 'ParametricCoefficients_20230721.mat' % output from R:patterns_of_neurodegeneration*.Rmd
        workdir
    end

    properties (Dependent)
        figdir
        label_permute
        N_bases_target
    end

    methods % GET
        function g = get.figdir(this)
            g = fullfile(this.workdir, ['baseline_', this.groups0{1}], 'results');
        end
        function g = get.label_permute(this)
            switch this.N_PATTERNS
                case 16
                    g = [1 2 13 16 17 18 19 20 21 22 3 4 5 6 7 8 9 10 11 12 14 15]; % ParamCoeff indices are 0, 1, 10, 11, 12, ...
                case 22
                case 24
                otherwise
                    error("mladni:ValueError", "%s: N_PATTERNS->%i", stackstr(), mladni.NMF.N_PATTERNS);
            end
        end
        function g = get.N_bases_target(this)
            g = this.N_PATTERNS;
        end
    end

    methods
        function [s1,s2] = call_apoe4(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_param));
            PC = ld.ParametricCoefficients20230721;

            apoe4 = PC(contains(PC.ParamCoefficients, "apoe4"), :);
            apoe4 = sortrows(apoe4, "ParamCoefficients");
            apoe4 = apoe4(this.label_permute, :);

            if opts.show
                disp(apoe4)
            end

            % Estimate
            P = apoe4.Estimate';
            SE = apoe4.StdErr';
            amin = min(P-SE, [], 'all');
            amax = max(P+SE, [], 'all');
            figure
            s1 = plot_with_stderr(this, P, SE, AxesMin=amin, AxesMax=amax, ...
                legend={'ApoE4'}, ti='Estimate \beta_{ApoE4}');
            saveFigures(this.figdir, closeFigure=true, prefix='E beta_apoe4');
    
            % PValue
            fprintf('NMFRadar.call_apoe4:\n')
            [~,~,~,P] = fdr_bh(apoe4.PValue', 0.05, 'dep', 'yes');
            P = asrow(P);
            amin = min(P, [], 'all');
            amax = max(P, [], 'all');
            if amin < 1e-15 || amax < 1e-15
                return
            end
            figure
            axes_scaling = repmat({'log'}, [this.N_PATTERNS 1]);
            s2 = plot(this, P, AxesMin=amin, AxesMax=amax, ...
                legend={'ApoE4'}, ti='FDR p-value \beta_{ApoE4}', AxesScaling=axes_scaling);
            saveFigures(this.figdir, closeFigure=true, prefix='FDR p-value beta_apoe4');
        end
        function [s1,s2] = call_groups(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            CC = ld.CohortCoefficients20230721;

            assert(length(this.groups) == length(this.groupLabels))
            s1 = cell(1, length(this.groups));
            s2 = cell(1, length(this.groups));
            for ig = 2:length(this.groups) % cn is reference cohort/category

                try
                    T = CC(contains(CC.ParamCoefficients, this.groups{ig}), :);
                    T = sortrows(T, "ParamCoefficients");
                    T = T(this.label_permute, :);
        
                    if opts.show
                        disp(T)
                    end
        
                    % Estimate
                    P = T.Estimate';
                    SE = T.StdErr';
                    amin = min(P-SE, [], 'all');
                    amax = max(P+SE, [], 'all');
                    figure
                    s1{ig} = plot_with_stderr(this, P, SE, AxesMin=amin, AxesMax=amax, ...
                        legend=this.groupLabels(ig), ...
                        ti="Estimate \beta_{" + this.groupLabels{ig} + "}");
                    saveFigures(this.figdir, closeFigure=true, prefix=sprintf('E beta_%s', this.groups{ig}));
        
                    % PValue
                    fprintf('NMFRadar.call_groups:\n')
                    [~,~,~,P] = fdr_bh(T.PValue, 0.05, 'dep', 'yes');
                    P = asrow(P);
                    amin = min(P, [], 'all');
                    amax = max(P, [], 'all');
                    if amin < 1e-15 || amax < 1e-15
                        continue
                    end
                    figure
                    axes_scaling = repmat({'log'}, [this.N_PATTERNS 1]);
                    s2{ig} = plot(this, P, AxesMin=amin, AxesMax=amax, ...
                        legend=this.groupLabels(ig), ...
                        ti="FDR p-value \beta_{" + this.groupLabels{ig} + "}", ...
                        AxesScaling=axes_scaling);
                    saveFigures(this.figdir, closeFigure=true, prefix=sprintf('FDR p-value beta_%s', this.groups{ig}));
                catch ME
                    handwarning(ME)
                end
            end
        end
        function [s1,s2] = call_sex(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_param));
            PC = ld.ParametricCoefficients20230721;

            male = PC(contains(PC.ParamCoefficients, "sexM"), :);
            male = sortrows(male, "ParamCoefficients");
            male = male(this.label_permute, :);

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
                legend={'male'}, ti='Estimate \beta_{sex}');
            saveFigures(this.figdir, closeFigure=true, prefix='E beta_sex');

            % PValue
            fprintf('NMFRadar.call_sex:\n')
            [~,~,~,P] = fdr_bh(male.PValue', 0.05, 'dep', 'yes');
            P = asrow(P);
            amin = min(P, [], 'all');
            amax = max(P, [], 'all');
            figure
            axes_scaling = repmat({'log'}, [this.N_PATTERNS; 1]);
            s2 = plot(this, P, AxesMin=amin, AxesMax=amax, ...
                legend={'male'}, ti='FDR p-value \beta_{sex}', AxesScaling=axes_scaling);
            saveFigures(this.figdir, closeFigure=true, prefix='FDR p-value beta_sex');
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
            plot(this, sigma, ...
                AxesMin = dipmin(sigma), AxesMax = dipmax(sigma), ...
                ti='S.D. Pattern-Weighted FDG')
            saveFigures(this.figdir, closeFigure=true, prefix='patt_weighted_fdg_sigma');

            figure
            plot(this, mu, ...
                AxesMin = dipmin(mu), AxesMax = dipmax(mu), ...
                ti='Mean Pattern-Weighted FDG')
            saveFigures(this.figdir, closeFigure=true, prefix='patt_weighted_fdg_mu');

            figure
            plot(this, snr, ...
                AxesMin = dipmin(snr), AxesMax = dipmax(snr), ...
                ti='\mu/\sigma Pattern-Weighted FDG')
            saveFigures(this.figdir, closeFigure=true, prefix='patt_weighted_fdg_snr');

            figure
            plot(this, cov, ...    
                AxesMin = dipmin(cov), AxesMax = dipmax(cov), ...
                ti='\sigma/\mu Pattern-Weighted FDG')
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
            end

            Nbases = size(P,2);
            
            % Delete variable in workspace if exists
            if exist('s', 'var')
                delete(s);
            end
            
            % Spider plot
            s = spider_plot_class(P);
            
            % Spider plot properties
            s.AxesLabels = ...
                {'ADNI P1', 'ADNI P2', 'ADNI P3', 'ADNI P4', ...
                'ADNI P5', 'ADNI P6', 'ADNI P7', 'ADNI P8', ...
                'ADNI P9', 'ADNI P10', 'ADNI P11', 'ADNI P12', ...
                'ADNI P13', 'ADNI P14', 'ADNI P15', 'ADNI P16', ...
                'ADNI P17', 'ADNI P18', 'ADNI P19', 'ADNI P20', ...
                'ADNI P21', 'ADNI P22'};
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
                opts.legend {mustBeText} = this.mergeDx
                opts.ti {mustBeTextScalar}
                opts.AxesMin double = []
                opts.AxesMax double = []
                opts.AxesInterval double = []
                opts.AxesScaling cell = {}
            end

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
            s.AxesLabels = ...
                {'ADNI P1', 'ADNI P2', 'ADNI P3', 'ADNI P4', ...
                'ADNI P5', 'ADNI P6', 'ADNI P7', 'ADNI P8', ...
                'ADNI P9', 'ADNI P10', 'ADNI P11', 'ADNI P12', ...
                'ADNI P13', 'ADNI P14', 'ADNI P15', 'ADNI P16', ...
                'ADNI P17', 'ADNI P18', 'ADNI P19', 'ADNI P20', ...
                'ADNI P21', 'ADNI P22'};
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
            ld = load(this.matfile_cohort);
            CC = ld.CohortCoefficients20230721;
            NP = this.N_PATTERNS;
            h = figure;
            hold on
            
            plot(CC.Estimate(1:NP), CC.Estimate(23:44), LineStyle="none", Marker="o", MarkerSize=15)
            plot(CC.Estimate(1:NP), CC.Estimate(43:66), LineStyle="none", Marker="s", MarkerSize=15)
            plot(CC.Estimate(1:NP), CC.Estimate(45:66), LineStyle="none", Marker="s", MarkerSize=15)
            plot(CC.Estimate(1:NP), CC.Estimate(89:110), LineStyle="none", Marker="*", MarkerSize=15)

            switch mladni.NMF.N_PATTERNS
                case 16
                    indices = {1 2 11 12 13 14 15 16 17 18 19 20 3 21 22 4 5 6 7 8 9 10};
                otherwise
                    error("mladni:ValueError", "%s: N_PATTERNS->%i", stackstr(), mladni.NMF.N_PATTERNS);
            end
            labels = cellfun(@(x) sprintf('P%i', x), indices, UniformOutput=false);
            labelpoints(CC.Estimate(1:NP), CC.Estimate(23:44), labels, 'SE', 0.2, 1)
            labelpoints(CC.Estimate(1:NP), CC.Estimate(45:66), labels, 'SE', 0.2, 1)
            labelpoints(CC.Estimate(1:NP), CC.Estimate(89:110), labels, 'SE', 0.2, 1)
        end
        
        function this = NMFRadar(varargin)
            this.workdir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
