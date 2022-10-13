classdef NMFRegression2 
    %% Accomodating Aris' preferences for regressions from 2022/09/29.
    %  GAMS implemented with R package mgcv RRID:SCR_001905 as described in
    %  Pehlivanova, M., Wolf, D.H., Sotiras, A., Kaczkurkin, A.N., Moore, T.M., Ciric, R., Cook, P.A., 
    %  Garza, A.G. de L., Rosen, A.F.G., Ruparel, K., Sharma, A., Shinohara, R.T., Roalf, D.R., Gur, R.C., 
    %  Davatzikos, C., Gur, R.E., Kable, J.W., Satterthwaite, T.D., 
    %  2018. Diminished Cortical Thickness Is Associated with Impulsive Choice in Adolescence. 
    %  J. Neurosci. 38, 2471–2481. https://doi.org/10.1523/JNEUROSCI.2200-17.2018
    %  
    %  Created 29-Sep-2022 01:12:50 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.13.0.2049777 (R2022b) for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        selected_covariates
        table_covariates
        table_selected
        table_selected_mat
        table_version
        workdir
    end

    methods
        function build_Aris_tables(this)
            g = this.table_selected;
            g = g(~strcmp(g.CDGLOBAL, '-1'), :);
            g.CDGLOBAL(g.AmyloidStatus == 1) = strcat(g.CDGLOBAL(g.AmyloidStatus == 1), ' amy+');
            g.CDGLOBAL(g.AmyloidStatus == 0) = strcat(g.CDGLOBAL(g.AmyloidStatus == 0), ' amy-');
            g.CDGLOBAL(isnan(g.AmyloidStatus)) = strcat(g.CDGLOBAL(isnan(g.AmyloidStatus)), ' amy n/a');
            c = g.Components;
            for idx = 1:28
                g.Properties.VariableNames = {'AmyloidStatus' 'MergeAge' 'MergePtGender' 'MergeCdrsbBl' 'CDGLOBAL' sprintf('Component%i', idx)};
                g.(sprintf('Component%i', idx)) = c(:,idx);
                writetable(g, sprintf('table_selected_%s_comp%i.csv', this.table_version, idx));
            end
        end

        function this = NMFRegression2()
            this.workdir = '/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline4/NumBases28/components';
            cd(this.workdir);
            ld = load('mladni_FDGDemographics_table_covariates.mat');
            this.table_covariates = ld.t;
            this.selected_covariates = {'MergePtGender', 'MergeAge', 'MergeCdrsbBl', 'CDGLOBAL', 'AmyloidStatus', 'Components'};

            g = table;
            g_varnames = {};            
            for varname = this.table_covariates.Properties.VariableNames
                if contains(varname{1}, this.selected_covariates)
                    g = addvars(g, this.table_covariates.(varname{1}));
                    g_varnames = [g_varnames varname{1}]; %#ok<AGROW>
                end
            end
            g.Properties.VariableNames = g_varnames;
            g.CDGLOBAL = convertStringsToChars(string(g.CDGLOBAL)); % numeric col -> char col
            this.table_selected = g;

            this.table_version = '20220930';
            this.table_selected_mat = sprintf('table_selected_%s.mat', this.table_version);
            if ~isfile(this.table_selected_mat)
                save(this.table_selected_mat, 'g');
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
