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
        table_selected
        workdir
    end

    methods
        function build_Aris_tables(this)
            g = this.table_selected;
            g = g(~strcmp(g.CDGLOBAL, '-1'), :);
            c = g.Components;
            for idx = 1:28
                g.Properties.VariableNames = {'MergeAge' 'MergePtGender' 'AmyloidStatus' 'CDGLOBAL' 'CdrsbBl' sprintf('Component%i', idx)};
                g.(sprintf('Component%i', idx)) = c(:,idx);
                writetable(g, sprintf('table_selected_20220929_comp%i.csv', idx));
            end
        end

        function this = NMFRegression2()            
            this.workdir = '/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline4/NumBases28/components';
            cd(this.workdir);
            ld = load('table_selected_20220929.mat');
            this.table_selected = ld.g;
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
