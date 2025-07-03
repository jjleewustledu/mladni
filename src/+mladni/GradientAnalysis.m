classdef GradientAnalysis < handle
    %% line1
    %  line2
    %  
    %  Created 27-Sep-2024 16:03:44 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2689473 (R2024a) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        homedir = fullfile(getenv("HOME"), "PycharmProjects", "gradient_analysis")
    end

    methods
        function this = GradientAnalysis(varargin)
        end

        function ic = build_super_volume(this)
            masksdir = fullfile(this.homedir, "gradient_data", "masks");
            g = mglob(fullfile(masksdir, "volume_*.nii.gz"));
            g = natsort(g);
            ifc = mlfourd.ImagingFormatContext2(g(1));
            for gidx = 2:length(g)
                ifc_ = mlfourd.ImagingFormatContext2(g(gidx));
                ifc.img(:,:,:,gidx) = ifc_.img;
            end
            ifc.fileprefix = "volume_super";
            ic = mlfourd.ImagingContext2(ifc);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
