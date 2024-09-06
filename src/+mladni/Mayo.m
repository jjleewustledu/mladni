classdef Mayo < DataCuration & handle
    %% line1
    %  line2
    %  
    %  Created 16-Apr-2024 22:41:50 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2537033 (R2024a) for MACA64.  Copyright 2024 John J. Lee.
    
    properties (Dependent)
        data_home
        nmf_fdg_home
    end

    methods % GET
        function g = get.data_home(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "MAYO");
        end
        function g = get.nmf_fdg_home(this)
            g = fullfile(this.data_home, "NMF_FDG");
        end
    end

    methods
        function this = Mayo()
            this.demogr = mladni.AdniDemographics();
            this.nmf = mladni.NMF( ...
                data_home=this.data_home, ...
                N_patterns=this.N_patterns);
            this.nmfc = mladni.NMFCovariates();
            this.nmfh = mladni.NMFHierarchies(data_home=this.data_home);
            this.nmfr = mladni.NMFRadar();
        end

        function call(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
