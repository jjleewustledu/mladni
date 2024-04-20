classdef Plot < handle
    %% line1
    %  line2
    %  
    %  Created 20-Apr-2024 13:02:50 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2568132 (R2024a) Update 1 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        fontsize_lab = 20
        fontsize_axes = 20
    end

    methods
        function this = Plot(varargin)
        end

        function h = plotxy(this, x, y, opts)
            %% provides a consistent, basic, x-y plotting for publications.
            %
            %  Also saves {.fig,.png,.svg}.
            %  For plot h, h.Position := [1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy].
            %
            %  Args:
            %    this mladni.Plot
            %    x double
            %    y double
            %    opts.xlab {mustBeText} = ""
            %    opts.ylab {mustBeText} = ""
            %    opts.fileprefix {mustBeTextScalar} = stackstr(3)
            %    opts.fracx double = 0.5
            %    opts.Npx double = 3400
            %    opts.fracy double = 0.5
            %    opts.Npy double = 1440

            arguments
                this mladni.Plot
                x double
                y double
                opts.xlab {mustBeText} = ""
                opts.ylab {mustBeText} = ""
                opts.fileprefix {mustBeTextScalar} = ""
                opts.fracx double = 0.5
                opts.Npx double = 2330
                opts.fracy double = 0.5
                opts.Npy double = 1440
            end

            h = figure;
            plot( ...
                ascol(x), ascol(y), ...
                '-o', LineWidth=4, Color=[0.5,0.5,0.5], ...
                MarkerSize=18, MarkerEdgeColor=[0,0,0], MarkerFaceColor=[0.5,0.5,0.5]);
            %xlim([x(1), x(end)]);
            xlabel(opts.xlab, fontsize=this.fontsize_lab);
            ylabel(opts.ylab, fontsize=this.fontsize_lab);
            set(gca, fontsize=this.fontsize_axes);
            set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy]);
            if ~isemptytext(opts.fileprefix)
                saveas(gcf, opts.fileprefix+".fig");
                saveas(gcf, opts.fileprefix+".png");
                saveas(gcf, opts.fileprefix+".svg");
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
