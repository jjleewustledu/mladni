classdef NMFRegression2 < handle
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
        Ncomp = 28
        selected_covariates
        table_covariates
        table_selected
        table_selected_mat
        table_version
        workdir
    end

    methods

        function this = NMFRegression2()
            this.workdir = '/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline4/NumBases28/components';
            cd(this.workdir);
            ld = load('mladni_FDGDemographics_table_covariates.mat');
            this.table_covariates = ld.t;
            this.table_covariates = this.table_covariates(this.table_covariates.CDGLOBAL ~= -1, :);
            writetable(this.table_covariates, 'mladni_FDGDemographics_table_covariates.csv');
            this.selected_covariates = {'MergePtGender', 'MergeAge', 'MergeCdrsbBl', 'CDGLOBAL', 'AmyloidStatus', 'Components'};

            g = table;
            g_varnames = {};
            selected = this.selected_covariates;
            for idx = 1:length(selected)
                if contains(selected{idx}, this.table_covariates.Properties.VariableNames)
                    g = addvars(g, this.table_covariates.(selected{idx}));
                    g_varnames = [g_varnames selected{idx}]; %#ok<AGROW>
                end
            end
            g.Properties.VariableNames = g_varnames;
            if any(contains(this.selected_covariates, 'AmyloidStatus'))
                g.AmyloidStatus(0 == g.AmyloidStatus) = -1;
                g.AmyloidStatus(isnan(g.AmyloidStatus)) = 0;
            end
            if any(contains(this.selected_covariates, 'CDGLOBAL'))
                g.CDGLOBAL = convertStringsToChars(string(g.CDGLOBAL)); % numeric col -> char col
            end
            this.table_selected = g;

            this.table_version = char(datetime('now', Format='yyyyMMdd'));
            this.table_selected_mat = sprintf('table_selected_%s.mat', this.table_version);
            deleteExisting(this.table_selected_mat)
            save(this.table_selected_mat, 'g');
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
