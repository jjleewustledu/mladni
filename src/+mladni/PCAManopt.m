classdef PCAManopt < handle
    %% Computes a version of NMF (nonnegative matrix factorization) on data.
    %
    % obj = mladni.PCAManopt(X, r);
    % call(obj)
    %
    % Given a matrix X of size m by n, such that each column represents a
    % point in R^m, this computes W: an orthonormal basis of size m by r such
    % that the column space of W captures the points X as well as possible.
    % More precisely, the function attempts to compute W as the minimizer
    % over the Grassmann manifold (the set of linear subspaces) of:
    %
    %  f(W) = (1/n) Sum_{i = 1:n} dist(X(:, i), the space spanned by W)
    %       = (1/n) Sum_{i = 1:n} || W*W'*X(:, i) - X(:, i) ||
    %
    % The output cost represents the average distance achieved with the
    % returned W. Notice that norms are not squared, for robustness.
    %
    % In practice, because this function is nonsmooth, it is smoothed with a
    % pseudo-Huber loss function of parameter epsilon (noted e for short), and
    % the smoothing parameter is iteratively reduced (with warm starts):
    %
    %   f_e(W) = (1/n) Sum_{i = 1:n} l_e(|| W*W'*X(:, i) - X(:, i) ||)
    %
    %   with l_e(x) = sqrt(x^2 + e^2) - e (for e = 0, this is absolute value).
    %
    % The intermediate optimization of the smooth cost over the Grassmann
    % manifold is performed using the Manopt toolbox.
    %
    % Ideally, the non-outlier data should be centered. If not, this
    % pre-processing centers all the data, but bear in mind that outliers will
    % shift the center of mass too.
    % X = X - repmat(mean(X, 2), [1, size(X, 2)]);
    %
    % There are no guarantees that this code will return the optimal W.
    % This code is distributed to illustrate one possible way of optimizing
    % a nonsmooth cost function over a manifold, using Manopt with smoothing.
    % For practical use, the constants in the code would need to be tuned.
    %
    % This file is part of Manopt and is copyrighted. See the license file.
    %
    % Main author: Nicolas Boumal and Teng Zhang, May 2, 2014
    % Contributors:
    %
    % Change log:
    %
    %   March 4, 2015 (NB):
    %       Uses a pseudo-Huber loss rather than a Huber loss: this has the
    %       nice advantage of being smooth and simpler to code (no if's).
    %
    %   April 8, 2015 (NB):
    %       Built-in test data for quick tests; added comment about centering.
    %
    %   Aug  20, 2021 (XJ):
    %       Added AD to compute the egrad
    %  
    %   Created 13-May-2024 13:55:53 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %   Developed on Matlab 24.1.0.2578822 (R2024a) Update 2 for MACA64.  Copyright 2024 John J. Lee.

    properties
        do_plot
        verbose
    end
    
    properties (SetAccess = protected)
        X
        m  % # voxels
        r  % rank of W
        n  % # imaging samples

        W  % X \approx \hat{X} = W * W' * X
        n_iterations
        epsilon
        reduction % of epsilon
        cost
        info
        problem
        residual
    end

    methods
        function this = PCAManopt(X, r, opts)
            arguments
                X {mustBeNumeric} = []
                r {mustBeScalarOrEmpty} = 1
                opts.n_iterations {mustBeScalarOrEmpty} = 6
                opts.reduction {mustBeScalarOrEmpty} = 0.5;
                opts.do_plot logical = true
                opts.verbose logical = true
            end
            if isempty(X)
                X = this.randn_X();
            end

            assert(contains(which("grassmannfactory"), "grassmannfactory"))
            this.X = X;
            this.r = r;
            this.n_iterations = opts.n_iterations;
            this.reduction = opts.reduction;

            this.do_plot = opts.do_plot;
            this.verbose = opts.verbose;
        end

        function call(this)
            [this.W, this.cost] = this.robust_pca(this.X);
            this.residual = norm(this.X - this.W * this.W' * this.X) / norm(this.X);
            fprintf("W = \n")
            disp(this.W)
            fprintf("mean(X,2) = \n")
            disp(mean(this.X, 2))
            fprintf("cost = %g\n", this.cost)
            fprintf("residual = %g\n", this.residual)
            
            plot_checks(this)
            if this.do_plot
                plot(this)
            end
        end

        function plot(this)
            X_ = this.X;
            W_ = this.W;
            m_ = size(X_, 1);
            assert(m_ > 1)            

            % Compare to a standard PCA
            [Wpca, ~, ~] = svds(X_, this.r);

            if 2 == m_
                assert(1 == this.r)                
                figure;
                scatter(X_(1,:), X_(2,:));
                range = 1.1*[min(X_(1,:)), max(X_(1,:)); ...
                             min(X_(2,:)), max(X_(2,:))];
                hold on;
                plot(range(1,:), range(2,:) * W_(2) / W_(1), 'r');
                hold off;
                hold on;
                plot(range(1,:), range(2,:) * Wpca(2) / Wpca(1), 'k');
                hold off;
                xlabel("\mu = " + 1);
                ylabel("\mu = " + 2);
                legend('data points', 'Robust PCA fit', 'Standard PCA fit');
                return
            end

            for rho = 1:this.r
                figure
                tiledlayout(m_-1, m_-1)
                for y = 2:m_
                    for x = 1:(m_-1)
                        if x < y
                            nexttile(x + (y-2)*(m_-1));
                            scatter(X_(x,:), X_(y,:));
                            range = 1.1*[min(X_(x,:)), max(X_(x,:)); ...
                                         min(X_(y,:)), max(X_(y,:))];
                            hold on;
                            plot(range(1,:), range(2,:) * W_(y, rho) / W_(x, rho), 'r');
                            hold off;
                            hold on;
                            plot(range(1,:), range(2,:) * Wpca(y, rho) / Wpca(x, rho), 'k');
                            hold off;
                            xlabel("\mu = " + x);
                            ylabel("\mu = " + y);
                            legend("data points, \rho = " + rho, "Robust PCA fit", "Standard PCA fit");
                        end
                    end
                end
            end
        end

        function plot_checks(this)

            % Numerically check gradient and Hessian consistency.
            figure;
            checkgradient(this.problem);
            if isfield(this.problem, "ehess")
                figure;
                checkhessian(this.problem);
            end

            % Display some statistics.
            if ~isempty(this.info)
                figure;
                semilogy([this.info.iter], [this.info.gradnorm], '.-');
                xlabel('Iteration #');
                ylabel('Gradient norm');
                title('Convergence of the trust-regions algorithm on the manifold');
            end
        end
    end

    methods (Static)
        function X_ = randn_X()
            %% generate random data for illustration purposes.

            % Generate some data points aligned on a subspace
            X_ = rand(2, 1)*(1:30) + .05*randn(2, 30).*[(1:30);(1:30)];
            % And add some random outliers to the mix
            P = randperm(size(X_, 2));
            outliers = 10;
            X_(:, P(1:outliers)) = 30*randn(2, outliers);
            % Center the data
            % X = X - repmat(mean(X, 2), [1, size(X, 2)]);
        end

        function X_ = rand_X(opts)
            %% generate nonnegative random data for illustration purposes.

            arguments
                opts.m {mustBeScalarOrEmpty} = 2
                opts.n {mustBeScalarOrEmpty} = 30
                opts.outliers {mustBeScalarOrEmpty} = 10
            end

            % Generate some data points aligned on a subspace
            X_ = rand(opts.m, 1)*(1:opts.n) + .05*rand(opts.m, opts.n).*repmat(1:opts.n, [opts.m, 1]);
            % And add some random outliers to the mix
            P = randperm(size(X_, 2));
            X_(:, P(1:opts.outliers)) = opts.n*rand(opts.m, opts.outliers);
        end

        function self_test()

            % Verify that Manopt was indeed added to the Matlab path.
            if isempty(which('spherefactory'))
                error(['You should first add Manopt to the Matlab path.\n' ...
     		       'Please run importmanopt.']);
            end

            % Generate the problem data.
            n_ = 1000;
            A = randn(n_);
            A = .5*(A+A');

            % Create the problem structure.
            manifold = spherefactory(n_);
            problem.M = manifold;

            % Define the problem cost function and its gradient.
            problem.cost  = @(x) -x'*(A*x);
            problem.egrad = @(x) -2*A*x;
            problem.ehess = @(x, xdot) -2*A*xdot;

            % Numerically check gradient and Hessian consistency.
            figure;
            checkgradient(problem);
            figure;
            checkhessian(problem);

            % Solve.
            [x, xcost, info] = trustregions(problem);  %#ok<ASGLU>

            % Display some statistics.
            figure;
            semilogy([info.iter], [info.gradnorm], '.-');
            xlabel('Iteration #');
            ylabel('Gradient norm');
            title('Convergence of the trust-regions algorithm on the sphere');
        end
    end

    %% PRIVATE

    methods (Access = private)
        function [W, cost] = robust_pca(this, X)

            % Prepare a Manopt problem structure for optimization of the given
            % cost (defined below) over the Grassmann manifold.
            [this.m, this.n] = size(X);
            manifold = grassmannfactory(this.m, this.r);
            problem_.M = manifold;
            problem_.cost = @robustpca_cost;
            % problem_.egrad = @robustpca_gradient;

            % An alternative way to compute the egrad is to use automatic
            % differentiation provided in the deep learning toolbox (slower).
            % Call manoptAD to automatically obtain the egrad
            problem_ = manoptAD(problem_, 'hess');

        	% Do classical PCA for the initial guess.
        	% This is just one idea: it is not necessarily useful or ideal.
            % Using a random initial guess, and starting over for a few different
            % ones is probably much better. For this example, we keep it simple.
            [W, ~, ~] = svds(X, this.r);


        	% Iteratively reduce the smoothing constant epsilon and optimize
        	% the cost function over Grassmann.
            this.epsilon = 1;
        	options.verbosity = 3; % Change this number for more or less output
            options.maxiter = 1e4;
            warning('off', 'manopt:getHessian:approx');

            for iter = 1 : (this.n_iterations - 1)
                W = trustregions(problem_, W, options);
                this.epsilon = this.epsilon * this.reduction;
            end
            [W, ~, info_] = trustregions(problem_, W, options);
            warning('on', 'manopt:getHessian:approx');

        	% Return the cost as the actual sum of distances, not smoothed.
        	this.epsilon = 0;
        	cost = robustpca_cost(W);

            % Store problem & info
            this.problem = problem_;
            this.info = info_;



            %% INNER FUNCTIONS

            % Smoothed cost
            function value = robustpca_cost(W)

                vecs = W*(W'*this.X) - this.X;
                sqnrms = sum(vecs.^2, 1);
                vals = sqrt(sqnrms + this.epsilon^2) - this.epsilon;
                value = mean(vals);

            end

            % Euclidean gradient of the smoothed cost (it will be transformed into
            % the Riemannian gradient automatically by Manopt).
            function G = robustpca_gradient(W)

        		% Note that the computation of vecs and sqnrms is redundant
        		% with their computation in the cost function. To speed
        		% up the code, it would be wise to use the caching capabilities
        		% of Manopt (the store structure). See online documentation.
        		% It is not done here to keep the code a bit simpler.

                % W = abs(W);  % produces negative curvature

                WtX = W'*this.X;
                vecs = W*WtX - this.X;  % sign reversal produce negative curvature
                sqnrms = sum(vecs.^2, 1);                

                % This explicit loop is a bit slow: the code below is equivalent
                % and faster to compute the gradient.
                % G = zeros(m, r);
                % for i=1:n
                %     G = G + (1/sqrt(sqnrms(i) + epsilon^2)) * vecs(:,i) * WtX(:,i)';
                % end
                % G = G/n;
                G = mean( ...
                    multiscale(1./sqrt(sqnrms + this.epsilon^2), ...
                    multiprod(reshape(vecs, [this.m, 1, this.n]), ...
                    multitransp(reshape(WtX, [this.r, 1, this.n])))), 3);
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
