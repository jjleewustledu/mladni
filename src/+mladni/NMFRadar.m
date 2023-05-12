classdef NMFRadar
    %% line1
    %  line2
    %  
    %  Created 16-Feb-2023 22:57:16 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        groups = { ...
            'Intercept', 'preclinical', 'emci', 'lmci', 'cdr_gt_0p5_apos'}
        groupLabels = {...
            'Intercept', 'Pre-clinical', 'Early MCI' 'Late MCI', 'Dementia'}
        groups0 = { ...
            'cn', 'preclinical', ...
            'cdr_0p5_apos_emci', 'cdr_0p5_apos_lmci', ...
            'cdr_gt_0p5_apos'}
        mergeDx = { ...
            'CN' 'pre-clinical' 'early MCI' 'late MCI' 'AD'}

        bases0 = [16 30 34 32 38]
        label_permute = [1,2,9,10,11,12,13,14,15,16,3,4,5,6,7,8]
        matfile0 = 'mladni_FDGDemographics_table_covariates_on_cn.mat'
        matfile_cohort = 'CohortCoefficients.mat' % output from R:patterns_of_neurodegeneration*.Rmd
        matfile_param = 'ParametricCoefficients.mat' % output from R:patterns_of_neurodegeneration*.Rmd
        workdir
    end

    properties (Dependent)
        N_bases_target
    end

    methods % GET
        function g = get.N_bases_target(this)
            g = this.bases0(1);
        end
    end

    methods
        function [s1,s2] = call_apoe4(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_param));
            PC = ld.ParametricCoefficients;

            apoe4 = PC(contains(PC.ParametricCoefficients, "apoe4"), :);
            apoe4 = sortrows(apoe4, "ParametricCoefficients");
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
            saveFigures(closeFigure=true, prefix='E beta_apoe4', ext='.svg')
    
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
            axes_scaling = repmat({'log'}, [16 1]);
            s2 = plot(this, P, AxesMin=amin, AxesMax=amax, ...
                legend={'ApoE4'}, ti='FDR p-value \beta_{ApoE4}', AxesScaling=axes_scaling);
            saveFigures(closeFigure=true, prefix='FDR p-value beta_apoe4', ext='.svg')
        end
        function [s1,s2] = call_groups(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_cohort));
            CC = ld.CohortCoefficients;

            assert(length(this.groups) == length(this.groupLabels))
            s1 = cell(1, length(this.groups));
            s2 = cell(1, length(this.groups));
            for ig = 1:length(this.groups) % cn is reference cohort/category

                T = CC(contains(CC.ParametricCoefficients, this.groups{ig}), :);
                T = sortrows(T, "ParametricCoefficients");
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
                saveFigures(closeFigure=true, prefix=sprintf('E beta_%s', this.groups{ig}), ext='.svg')
    
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
                axes_scaling = repmat({'log'}, [16 1]);
                s2{ig} = plot(this, P, AxesMin=amin, AxesMax=amax, ...
                    legend=this.groupLabels(ig), ...
                    ti="FDR p-value \beta_{" + this.groupLabels{ig} + "}", ...
                    AxesScaling=axes_scaling);
                saveFigures(closeFigure=true, prefix=sprintf('FDR p-value beta_%s', this.groups{ig}), ext='.svg')
            end
        end
        function [s1,s2] = call_sex(this, opts)
            arguments
                this mladni.NMFRadar
                opts.show logical = false
            end

            ld = load(fullfile(this.workdir, this.matfile_param));
            PC = ld.ParametricCoefficients;

            male = PC(contains(PC.ParametricCoefficients, "sexM"), :);
            male = sortrows(male, "ParametricCoefficients");
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
            saveFigures(closeFigure=true, prefix='E beta_sex', ext='.svg')

            % PValue
            fprintf('NMFRadar.call_sex:\n')
            [~,~,~,P] = fdr_bh(male.PValue', 0.05, 'dep', 'yes');
            P = asrow(P);
            amin = min(P, [], 'all');
            amax = max(P, [], 'all');
            figure
            axes_scaling = repmat({'log'}, [16 1]);
            s2 = plot(this, P, AxesMin=amin, AxesMax=amax, ...
                legend={'male'}, ti='FDR p-value \beta_{sex}', AxesScaling=axes_scaling);
            saveFigures(closeFigure=true, prefix='FDR p-value beta_sex', ext='.svg')
        end        
        function this = call0(this)
            mu = nan(length(this.mergeDx), this.N_bases_target);
            sigma = nan(length(this.mergeDx), this.N_bases_target);
            snr = nan(length(this.mergeDx), this.N_bases_target);
            cov = nan(length(this.mergeDx), this.N_bases_target);
            grad = nan(length(this.mergeDx), this.N_bases_target);
            for c = 1:length(this.mergeDx)
                ld = load(fullfile(this.workdir, ...
                    sprintf('baseline_%s', this.groups0{c}), ...
                    sprintf('NumBases%i', this.N_bases_target), ...
                    'components', this.matfile0));

                mu(c,:) = mean(ld.t.Components, 1);
                sigma(c,:) = std(ld.t.Components, 1);
                snr(c,:) = mu(c,:)./sigma(c,:);
                cov(c,:) = sigma(c,:)./mu(c,:);
                grad(c,:) = (mu(c,:) - mu(1,:))./mu(1,:);
            end

            figure
            plot(this, grad, ...
                AxesMin = -0.2, AxesMax = 0.025, ...
                ti='\Delta_{CN}/\mu_{CN} Pattern Coefficients for Diagnostic Cohorts')


            return


            figure
            plot(this, mu, ...
                AxesMin = 0.7, AxesMax = 1.35, ...
                ti='Mean Pattern Coefficients for Diagnostic Cohorts')

            figure
            plot(this, sigma, ...
                AxesMin = 0.03, AxesMax = 0.17, ...
                ti='S.D. Pattern Coefficients for Diagnostic Cohorts')

            figure
            plot(this, snr, ...
                AxesMin = 4, AxesMax = 25, ...
                ti='\mu/\sigma Pattern Coefficients for Diagnostic Cohorts')

            figure
            plot(this, cov, ...
                AxesMin = 0.03, AxesMax = 0.16, ...
                ti='\sigma/\mu Pattern Coefficients for Diagnostic Cohorts')
        end
        function s = plot(this, P, opts)
            %% spider_plot_class_examples.m Example 5: Excel-like radar charts.

            arguments
                this mladni.NMFRadar %#ok<INUSA> 
                P double
                opts.legend {mustBeText} = this.mergeDx
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
                'ADNI P13', 'ADNI P14', 'ADNI P15', 'ADNI P16'};
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
            s.AxesFontSize = 13;
            s.LabelFontSize = 14;
            s.AxesColor = [0.8, 0.8, 0.8];
            s.AxesLabelsEdge = 'none';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            s.LegendLabels = opts.legend;
            s.LegendHandle.Location = 'northeastoutside';
            s.LegendHandle.FontSize = 13;
            title(opts.ti, 'FontSize', 16);
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
            s.AxesLabels = { ...
                'ADNI P1', 'ADNI P2', 'ADNI P3', 'ADNI P4', ...
                'ADNI P5', 'ADNI P6', 'ADNI P7', 'ADNI P8', ...
                'ADNI P9', 'ADNI P10', 'ADNI P11', 'ADNI P12', ...
                'ADNI P13', 'ADNI P14', 'ADNI P15', 'ADNI P16'};
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
            s.AxesFontSize = 13;
            s.LabelFontSize = 14;
            s.AxesColor = [0.8, 0.8, 0.8];
            s.AxesLabelsEdge = 'none';
            s.AxesRadial = 'off';
            if ~isempty(opts.AxesScaling)
                s.AxesScaling = opts.AxesScaling;
            end
            
            s.LegendLabels = opts.legend;
            s.LegendHandle.Location = 'northeastoutside';
            s.LegendHandle.FontSize = 13;
            title(opts.ti, 'FontSize', 16);
        end
        function this = NMFRadar(varargin)
            this.workdir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
