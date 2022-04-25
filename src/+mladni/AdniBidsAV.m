classdef AdniBidsAV < mladni.AdniBids
    %% 
    %  
    %  Created 07-Apr-2022 14:32:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    
    methods
        function this = AdniBidsAV(varargin)
            %% ADNIBIDSAV 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            this = this@mladni.AdniBids(varargin{:})
            
            ip = inputParser;
            addParameter(ip, "arg1", [], @(x) false)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
