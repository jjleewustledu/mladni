classdef Plot < handle
    %% line1
    %  line2
    %  
    %  Created 20-Apr-2024 13:02:50 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2568132 (R2024a) Update 1 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        fh_central_tendency = @median  % @mean
        fontsize_lab = 20
        fontsize_axes = 20
        fontsize_scale = 2.5
        fractional_grey = 0.33
        interp = @pchip
        line_width = 5
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
            %    opts.fileprefix {mustBeTextScalar} = ""
            %    opts.fracx double = 0.5
            %    opts.Npx double = 2330
            %    opts.fracy double = 0.5
            %    opts.Npy double = 1440
            %    opts.do_interp logical = true

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
                opts.do_interp logical = true
            end

            if opts.do_interp
                dx = (x(end) - x(1)) / 10000;
                x_mak = x(1):dx:x(end);
                y_mak = this.interp(x, y, x_mak);
            end

            h = plot( ...
                ascol(x_mak), ascol(y_mak), LineWidth=4, Color=[0.5,0.5,0.5]);
                % '-o', LineWidth=4, Color=[0.5,0.5,0.5], ...
                % MarkerSize=18, MarkerEdgeColor=[0,0,0], MarkerFaceColor=[0.5,0.5,0.5]);
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

        function [h,p] = rm_plot(this, span_model, reproducibility, opts)
            %% provides a consistent, basic, x-y plotting for publications.
            %
            %  Also saves {.fig,.png,.svg}.
            %  For plot h, h.Position := [1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy].

            arguments
                this mladni.Plot
                span_model double
                reproducibility double
                opts.colours {mustBeNumeric} = this.fractional_grey * [1, 1, 1]
                opts.do_interp logical = true
                opts.line_width double = this.line_width
                opts.line_colour {mustBeNumeric} = this.fractional_grey * [1, 1, 1]
                opts.do_flip logical = true
            end

            marker_size = 200;

            if opts.do_flip
                span_model = flip(span_model);

                % confidence bounds
                if ismatrix(reproducibility)  % N_reps x spans
                    iqrs = iqr(reproducibility, 1);
                    reproducibility = median(reproducibility, 1);

                    if opts.do_interp
                        dspan = (span_model(end) - span_model(1)) / 100000;
                        span_model_mak = span_model(1):dspan:span_model(end);
                        confidence_x = [span_model_mak, flip(span_model_mak)];
                        iqrs_mak = this.interp(span_model, iqrs, span_model_mak);
                        reproducibility_mak = this.interp(span_model, reproducibility, span_model_mak);
                        confidence_y = [reproducibility_mak - iqrs_mak, flip(reproducibility_mak + iqrs_mak)];
                    else
                        confidence_x = [span_model, flip(span_model)];
                        confidence_y = [reproducibility - iqrs, flip(reproducibility + iqrs)];
                    end

                    p = fill(confidence_y, confidence_x, opts.colours, ...
                        FaceColor=opts.colours, ...
                        EdgeColor="none", ...
                        FaceAlpha=0.15);
                end

                % interpolated trendline
                if opts.do_interp
                    dspan = (span_model(end) - span_model(1)) / 100000;
                    span_model_mak = span_model(1):dspan:span_model(end);
                    reproducibility_mak = this.interp(span_model, reproducibility, span_model_mak);
                    h = plot( ...
                        ascol(reproducibility_mak), ascol(span_model_mak), ...
                        LineWidth=opts.line_width, ...
                        Color=opts.line_colour);
                    return
                end

                % scatter plot
                h = scatter(reproducibility, span_model, ...
                    'MarkerFaceColor', opts.line_colour, ...
                    'MarkerEdgeColor', [0 0 0], ...
                    'MarkerFaceAlpha', 1, ...
                    'SizeData', marker_size, ...
                    'LineWidth', 2);

                return
            end

            % confidence bounds
            if ismatrix(reproducibility)  % N_reps x spans
                iqrs = iqr(reproducibility, 1);
                reproducibility = median(reproducibility, 1);

                if opts.do_interp
                    dspan = (span_model(end) - span_model(1)) / 100000;
                    span_model_mak = span_model(1):dspan:span_model(end);
                    confidence_x = [span_model_mak, flip(span_model_mak)];
                    iqrs_mak = this.interp(span_model, iqrs, span_model_mak);
                    reproducibility_mak = this.interp(span_model, reproducibility, span_model_mak);
                    confidence_y = [reproducibility_mak - iqrs_mak, flip(reproducibility_mak + iqrs_mak)];
                else
                    confidence_x = [span_model, flip(span_model)];
                    confidence_y = [reproducibility - iqrs, flip(reproducibility + iqrs)];
                end

                p = fill(confidence_x, confidence_y, opts.colours, ...
                    FaceColor=opts.colours, ...
                    EdgeColor="none", ...
                    FaceAlpha=0.15);
            end

            % interpolated trendline
            if opts.do_interp
                dspan = (span_model(end) - span_model(1)) / 100000;
                span_model_mak = span_model(1):dspan:span_model(end);
                reproducibility_mak = this.interp(span_model, reproducibility, span_model_mak);
                h = plot( ...
                    ascol(span_model_mak), ascol(reproducibility_mak), ...
                    LineWidth=opts.line_width, ...
                    Color=opts.line_colour);
                return
            end

            % scatter plot
            h = scatter(span_model, reproducibility, ...
                'MarkerFaceColor', opts.line_colour, ...
                'MarkerEdgeColor', [0 0 0], ...
                'MarkerFaceAlpha', 1, ...
                'SizeData', marker_size, ...
                'LineWidth', 2);
        end

        function h = rm_raincloud(this, D, opts)
            %% pattern_nums is numeric.
            %  pattern_names is text.
            %  D is data ~ N_reps x K_patterns.

            arguments
                this mladni.Plot
                D 
                opts.rain_colours = this.fractional_grey * [1, 1, 1]
                opts.pdf_colours = this.fractional_grey * [1, 1, 1]
                opts.do_interp logical = true
                opts.do_plot_dots logical = false
                opts.line_width double = this.line_width
                opts.line_colour = this.fractional_grey * [1, 1, 1]
            end

            %% aufbau cell-array data

            if iscell(D)
                N_series = length(D);
                N_spans = min(cellfun(@(x) size(x, 2), D));
                data = cell(N_spans, N_series);
                for span = 1:N_spans
                    for series = 1:N_series
                        try
                            data{span, series} = asrow(D{series}(:, span));  % aufbau rows of bootstrap repetitions
                        catch ME
                            disp(ME.message)
                        end
                    end
                end
            else
                N_spans = size(D, 2);
                data = cell(N_spans, 1);
                for span = 1:N_spans
                    data{span, 1} = asrow(D(:, span));
                end
            end

            %% build rainclouds

            h = this.rm_raincloud_(data, ...
                rain_colours=opts.rain_colours, ...
                pdf_colours=opts.pdf_colours, ...
                do_interp=opts.do_interp, ...
                do_plot_dots=opts.do_plot_dots, ...
                line_width=opts.line_width, ...
                line_colour=opts.line_colour);
            set(gca, 'YTickLabel', 2*N_spans:-2:2);
            set(gca, 'XLim', [0, 1.1]);
        end

        function h = rm_raincloud_(this, data, opts)
            % Use like: h = rm_raincloud(data, colours, plot_top_to_bottom, density_type, bandwidth)
            % Where 'data' is an M x N cell array of N data series and M measurements
            % And 'colours' is an N x 3 array defining the colour to plot each data series
            % plot_top_to_bottom: Default plots left-to-right, set to 1 to rotate.
            % density_type: 'ks' (default) or 'RASH'. 'ks' uses matlab's inbuilt 'ksdensity' to
            % determine the shape of the rainclouds. 'RASH' will use the 'rst_RASH'
            % method from Cyril Pernet's Robust Stats toolbox, if that function is on
            % your matlab path.
            % bandwidth: If density_type == 'ks', determines bandwidth of density estimate
            % h is a cell array containing handles of the various figure parts:
            % h.p{i,j} is the handle to the density plot from data{i,j}
            % h.s{i,j} is the handle to the 'raindrops' (individual datapoints) from data{i,j}
            % h.m(i,j) is the handle to the single, large dot that represents mean(data{i,j})
            % h.l(i,j) is the handle for the line connecting h.m(i,j) and h.m(i+1,j)

            %% TO-DO:
            % Patch can create colour gradients using the 'interp' option to 'FaceColor'. Allow this?

            arguments
                this mladni.Plot
                data cell
                opts.raindrop_size = 50;
                opts.rain_colours 
                opts.pdf_colours 
                opts.plot_top_to_bottom logical = false;  % left-to-right plotting by default
                opts.density_type {mustBeTextScalar} = 'ks';  % use 'ksdensity' to create cloud shapes
                opts.bandwidth {mustBeNumeric} = [];  % let the function specify the bandwidth
                opts.do_interp logical = true;
                opts.do_plot_dots logical = false;
                opts.line_width double = this.line_width
                opts.line_colour = this.fractional_grey * [1, 1, 1]
            end
            opts.rain_colours = ensureCell(opts.rain_colours);
            opts.pdf_colours = ensureCell(opts.pdf_colours);
            opts.line_colour = ensureCell(opts.line_colour);

            % make sure we have correct number of colours
            % assert(all(size(colours) == [n_series 3]), 'number of colors does not match number of plot series');

            %% check dimensions of data

            [n_plots_per_series, n_series] = size(data);
            if 1 == size(opts.rain_colours, 1)
                opts.rain_colours = repmat(opts.rain_colours, [n_plots_per_series, 1]);
            end
            if 1 == size(opts.pdf_colours, 1)
                opts.pdf_colours = repmat(opts.pdf_colours, [n_plots_per_series, 1]);
            end

            %% Calculate properties of density plots

            % Probably okay to hard-code this as it just determines the granularity of
            % the density estimate
            density_granularity = 200;

            n_bins = repmat(density_granularity, n_plots_per_series, n_series);

            % calculate kernel densities
            for i = 1:n_plots_per_series
                for j = 1:n_series
                    try
                        switch opts.density_type

                            case 'ks'

                                % compute density using 'ksdensity'
                                [ks{i, j}, x{i, j}] = ksdensity(data{i, j}, 'NumPoints', n_bins(i, j), 'bandwidth', opts.bandwidth);

                            case 'rash'

                                % check for rst_RASH function (from Robust stats toolbox) in path, fail if not found
                                assert(exist('rst_RASH', 'file') == 2, 'Could not compute density using RASH method. Do you have the Robust Stats toolbox on your path?');

                                % compute density using RASH
                                [x{i, j}, ks{i, j}] = rst_RASH(data{i, j});

                                % override default 'n_bins' as rst_RASH determines number of bins
                                n_bins(i, j) = size(ks{i, j}, 2);
                        end

                        % Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
                        q{i, j}     = (1:n_bins(i, j) - 1)';
                        faces{i, j} = [q{i, j}, q{i, j} + 1, q{i, j} + n_bins(i, j) + 1, q{i, j} + n_bins(i, j)];

                    catch ME
                        disp(ME.message)
                    end
                end
            end

            % determine spacing between plots
            nontrivial  = cellfun(@(x) ~isempty(x), ks);
            spacing     = 2 * mean(mean(cellfun(@max, ks(nontrivial))));
            ks_offsets  = [0:n_plots_per_series-1] .* spacing;

            % flip so first plot in series is plotted on the *top*
            ks_offsets  = fliplr(ks_offsets);

            % calculate patch vertices from kernel density
            for i = 1:n_plots_per_series
                for j = 1:n_series
                    try

                        verts{i, j} = [x{i, j}', ks{i, j}' + ks_offsets(i); x{i, j}', ones(n_bins(i, j), 1) * ks_offsets(i)];
                        verts{i, j} = [x{i, j}', ks{i, j}' + ks_offsets(i); x{i, j}', ones(n_bins(i, j), 1) * ks_offsets(i)];

                    catch ME
                        disp(ME.message)
                    end
                end
            end


            %% jitter for the raindrops

            jit_width = spacing / 8;

            for i = 1:n_plots_per_series
                for j = 1:n_series
                    try
                       
                        jit{i,j} = jit_width + rand(1, length(data{i,j})) * jit_width;

                    catch ME
                        disp(ME.message)
                    end
                end
            end

            %% means/medians (for mean/median dots)

            cell_means = cellfun(this.fh_central_tendency, data);

            %% plot
            % note - we *could* plot everything here in one big loop, but then
            % different figure parts would overlay each other in a silly way.

            hold on

            % PDFs and raindrops
            for i = 1:n_plots_per_series
                for j = 1:n_series
                    try

                        % patch PDFs
                        h.p{i, j} = patch('Faces', faces{i, j}, 'Vertices', verts{i, j}, 'FaceVertexCData', opts.pdf_colours{i, j}, 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 0.25);

                        % scatter raindrops
                        h.s{i, j} = scatter(data{i, j}, -jit{i, j} + ks_offsets(i), 'MarkerFaceColor', opts.rain_colours{i, j}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.25, 'SizeData', opts.raindrop_size);

                    catch ME
                        disp(ME.message)
                    end
                end
            end

            % trendlines
            if opts.do_interp
                for j = 1:n_series

                    d_ks_offsets = (ks_offsets(end) - ks_offsets(1)) / 100000;
                    ks_offsets_mak = ks_offsets(1):d_ks_offsets:ks_offsets(end);
                    cell_means_mak = this.interp(ks_offsets, cell_means(:,j)', ks_offsets_mak);
                    h.interp = line(cell_means_mak, ks_offsets_mak, 'LineWidth', opts.line_width, 'Color', opts.line_colour{j});

                end
            else
                for i = 1:n_plots_per_series - 1 % We have n_plots_per_series-1 lines because lines connect pairs of points
                    for j = 1:n_series
                        try

                            h.l(i, j) = line(cell_means([i i+1], j), ks_offsets([i i+1]), 'LineWidth', 4, 'Color', opts.line_colour{j});

                        catch ME
                            disp(ME.message)
                        end
                    end
                end
            end

            % markers on trendlines
            if opts.do_plot_dots
                for i = 1:n_plots_per_series
                    for j = 1:n_series
                        try

                            h.m(i, j) = scatter(cell_means(i, j), ks_offsets(i), ...
                                'MarkerFaceColor', opts.pdf_colours(i, :), ...
                                'MarkerEdgeColor', [0 0 0], ...
                                'MarkerFaceAlpha', 1, ...
                                'SizeData', opts.raindrop_size * 2, ...
                                'LineWidth', 2);

                        catch ME
                            disp(ME.message)
                        end
                    end
                end
            end

            %% clear up axis labels

            % 'YTick', likes values that *increase* as you go up the Y-axis, but we plot the first
            % raincloud at the top. So flip the vector around
            set(gca, 'YTick', fliplr(ks_offsets));
            set(gca, 'YTickLabel', n_plots_per_series:-1:1);

            %% determine plot rotation
            % default option is left-to-right
            % opts.plot_top_to_bottom can be set to 1
            % NOTE: Because it's easier, we actually prepare everything plotted
            % top-to-bottom, then - by default - we rotate it here. That's why the
            % logical is constructed the way it is.

            % rotate and flip
            if ~opts.plot_top_to_bottom
                view([90 -90]);
                axis ij
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
