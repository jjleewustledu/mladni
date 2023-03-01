classdef Neurosynth
    %% line1
    %  line2
    %  
    %  Created 02-Feb-2023 23:20:54 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    methods (Static)
        function add_all()
            ga = glob(fullfile('topic-*', 'v4-topics-50_*association*.nii'))';
            gu = glob(fullfile('topic-*', 'v4-topics-50_*uniformity*.nii'))';
            assert(length(ga) == length(gu))
            N = length(ga);

            assoc = mlfourd.ImagingContext2(ga{1});
            unif = mlfourd.ImagingContext2(gu{1});
            for n = 2:N
                assoc = assoc + mlfourd.ImagingContext2(ga{n});
                assoc.fileprefix = 'v4-topics-50_association-test_z_FDR_0.01';
                unif = unif + mlfourd.ImagingContext2(gu{n});
                unif.fileprefix = 'v4-topics-50_uniformity-test_z_FDR_0.01';
            end
            assoc.save
            unif.save
        end
    end

    methods
        function this = Neurosynth(varargin)
            %% NEUROSYNTH 
            
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
