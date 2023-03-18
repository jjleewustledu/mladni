classdef NMFRadar
    %% line1
    %  line2
    %  
    %  Created 16-Feb-2023 22:57:16 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        groups = { ...
            'cn', 'preclinical', 'emci', 'lmci', 'cdr_gt_0p5_apos'}
        groups0 = { ...
            'cn', 'preclinical', ...
            'cdr_0p5_apos_emci', 'cdr_0p5_apos_lmci', ...
            'cdr_gt_0p5_apos'}
        mergeDx = { ...
            'CN' 'pre-clinical' 'early MCI' 'late MCI' 'AD'}

        bases0 = [16 30 34 32 38]
        matfile0 = 'mladni_FDGDemographics_table_covariates_on_cn.mat'
        matfile = 'CohortCoefficients.mat'
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
        function this = call(this)
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
        function h = plot(this, P, opts)
            %% spider_plot_class_examples.m Example 5: Excel-like radar charts.

            arguments
                this mladni.NMFRadar
                P double
                opts.ti {mustBeTextScalar}
                opts.AxesMin double
                opts.AxesMax double
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
            %s.AxesInterval = 5;
            s.AxesPrecision = 2;
            s.AxesDisplay = 'one';
            s.AxesLimits = [opts.AxesMin*ones(1,Nbases); opts.AxesMax*ones(1,Nbases)];
            s.FillOption = 'on';
            s.FillTransparency = 0.1;
            s.Color = [139, 0, 0; 240, 128, 128]/255;
            s.LineWidth = 2;
            s.Marker = 'none';
            s.AxesFontSize = 14;
            s.LabelFontSize = 10;
            s.AxesColor = [0.8, 0.8, 0.8];
            s.AxesLabelsEdge = 'none';
            s.AxesRadial = 'off';
            
            s.LegendLabels = this.mergeDx;
            title(opts.ti, 'FontSize', 14);

        end
        function this = NMFRadar(varargin)
            this.workdir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
