classdef NMFManopt < handle
    %% Computes a version of NMF (nonnegative matrix factorization) on data,
    %  using manopt's example truncated_svd() as an algorithmic source.
    %
    %  Usage:
    %  >> obj = mladni.NMFManopt(X, r);
    %  >> call(obj)
    %
    %  Args:
    %    X {mustBeNumeric}:  A real matrix of size mxn.
    %    r {mustBeScalar}:  An integer r <= min(m, n).
    %
    %  truncated_svd() is part of Manopt and is copyrighted. See the license file at manopt.org
    %
    %  Main author: Nicolas Boumal, July 5, 2013
    %
    %  Change log:
    %
    %  Aug. 23, 2021 (XJ):
    %      Added AD to compute the egrad and the ehess
    %  
    %  File ceated 13-May-2024 13:55:53 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2578822 (R2024a) Update 2 for MACA64.  Copyright 2024 John J. Lee.

    properties
        do_AD
        do_force_nonneg
        do_plot
        do_yang_oja
        solver_name
        verbose
    end
    
    properties (SetAccess = protected)
        X
        r  % rank of W
        
        batchsize
        cost_result
        epsilon0
        info
        manifold_factory  % handle to manopt factories
        maxiter
        n_iterations
        problem
        recon
        reduction  % of epsilon0
        result
        residual
        stepsize_lambda
        W0
    end

    properties (Dependent)
        m  % # voxels in row space of this.X
        n  % # imaging samples in col space of this.X

        W  % X \approx \hat{X} = W * W' * X
        H  % H = W' * X
    end

    methods
        function g = get.m(this)
            g = size(this.X, 1);
        end
        function g = get.n(this)
            g = size(this.X, 2);
        end
        function g = get.W(this)
            g = this.result;
            g(g < eps) = eps;
        end
        function g = get.H(this)
            g = this.W'*this.X;
        end
    end

    methods
        function this = NMFManopt(X, r, W0, manifold_factory, opts)
            arguments
                X double {mustBeNumeric} = []
                r double {mustBeScalarOrEmpty} = 1
                W0 double {mustBeNumeric} = []  % empty W0 will trigger random starting conditions for manopt
                manifold_factory function_handle = @grassmannfactory
                opts.batchsize {mustBeScalarOrEmpty} = floor(size(X,2)/10)
                opts.do_AD logical = true
                opts.do_force_nonneg logical = false
                opts.do_plot logical = true
                opts.do_yang_oja logical = false
                opts.epsilon0 double {mustBeScalarOrEmpty} = 0
                opts.maxiter {mustBeScalarOrEmpty} = 2e4
                opts.n_iterations {mustBeScalarOrEmpty} = 6
                opts.reduction double {mustBeScalarOrEmpty} = 0.5;
                opts.solver_name {mustBeTextScalar} = "stochasticgradient"
                opts.stepsize_lambda double {mustBeScalarOrEmpty} = 1
                opts.verbose logical = true
            end
            if isempty(X)
                [X, r] = this.rand_X();
            end
            if contains(opts.solver_name, "stochasticgradient")
                opts.do_AD = false;
            end

            assert(contains(which("grassmannfactory"), "grassmannfactory"))
            this.X = X;
            this.r = r;
            this.W0 = W0;
            this.manifold_factory = manifold_factory;

            this.batchsize = opts.batchsize;
            this.do_AD = opts.do_AD;
            this.do_force_nonneg = opts.do_force_nonneg;
            this.do_plot = opts.do_plot;
            this.do_yang_oja = opts.do_yang_oja;
            this.epsilon0 = opts.epsilon0;
            this.maxiter = opts.maxiter;
            this.n_iterations = opts.n_iterations;
            this.reduction = opts.reduction;
            this.solver_name = opts.solver_name;
            this.stepsize_lambda = opts.stepsize_lambda;
            this.verbose = opts.verbose;
        end

        function call(this)
            fprintf("================= start of %s =================\n", stackstr())
            this.nmf(this.X, this.r, this.do_AD, this.epsilon0);
            this.recon = this.W * this.H;
            this.residual = norm(this.recon - this.X, "fro") / norm(this.X, "fro");            
            fprintf("cost_value = %g\n", this.cost_result)
            fprintf("residual = %g\n", this.residual)
            fprintf("================= end of %s =================\n", stackstr())
            
            if this.do_plot
                plot(this)
            end
        end

        function call_svd(this)
            [this.W, S, V, this.info] = this.truncated_svd(this.X, this.r, this.do_AD);
            this.H = S * V';
            this.recon = this.W * this.H;
            this.residual = norm(this.X - this.W * this.H, "fro") / norm(this.X, "fro");            
            fprintf("cost_value = %g\n", this.cost_result)
            fprintf("residual = %g\n", this.residual)
            
            if this.do_plot
                plot(this)
            end
        end

        function plot(this)
            %%

            try
                % For our information, Manopt can also compute the spectrum of the
                % Riemannian Hessian on the tangent space at (any) X. Computing the
                % spectrum at the solution gives us some idea of the conditioning of
                % the problem. If we were to implement a preconditioner for the
                % Hessian, this would also inform us on its performance.
                %
                % Notice that if the optimization is performed on a product of Stiefel
                % manifolds instead of a product of Grassmannians, the double
                % invariance under the orthogonal group O(p) will appear as twice
                % p*(p-1)/2, thus p*(p-1) zero eigenvalues in the spectrum of the
                % Hessian. This means that the minimizers are not isolated, which
                % typically hinders convergence of second order algorithms.
                if this.problem.M.dim() < 1024
                    evs = hessianspectrum(this.problem, this.result);
                    stairs(sort(evs));
                    title(['Eigenvalues of the Hessian of the cost function ' ...
                        'at the solution']);
                end
            catch ME
                handwarning(ME)
            end

            % plots for stochasticgradient
            if isfield(this.info, "metric")
                figure

                % Plot the special metric recorded by options.statsfun
                semilogy([this.info.iter], [this.info.metric], '.-');
                xlabel('Iteration #');
                ylabel('Frobenius norm of WW^TX - X');
                title('Convergence of stochasticgradient');

                % Add to that plot a reference: the globally optimal value attained if
                % the true dominant singular vectors are computed.
                fprintf('Running svds... ');
                t = tic();
                [V, ~] = svds(this.X, this.r);
                fprintf('done: %g [s] (note: svd may be faster)\n', toc(t));
                hold all;
                bound = norm(this.X'*V, 'fro');
                plot([this.info.iter], bound*ones(size([this.info.iter])), '--');
                hold off;
                legend('Algorithm', 'SVD bound', 'Location', 'SouthEast');
            end

            % plot patterns
            figure
            hold on
            for rho = 1:this.r
                Wcol = smoothdata(this.W(:, rho), "movmean", 50);
                % Wcol = this.W(:, rho);
                if rho == 1
                    plot(asrow(Wcol), LineWidth=1);
                elseif rho == 2
                    plot(asrow(Wcol), LineWidth=0.5);
                else
                    plot(asrow(Wcol), ':');
                end
            end
            hold off
            xlabel("voxels")
            ylabel("suvr")
            title(sprintf("W(:, rho=1:%g)", this.r))
            legend("rho = " + (1:this.r))

            % plot recon
            figure 
            hold on
            plot(mean(this.X, 2), LineWidth=1)
            plot(mean(this.recon, 2), '--', LineWidth=2)
            hold off
            xlabel("voxels")
            ylabel("suvr")
            legend(["mean(X, 2)", "mean(recon, 2)"])
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
   
        function [W,H] = NNDSVD(A,k,flag)
            %% https://github.com/asotiras/brainparts/blob/master/NNDSVD.m
            % This function implements the NNDSVD algorithm described in [1] for
            % initialization of Nonnegative Matrix Factorization Algorithms.
            %
            % [W,H] = nndsvd(A,k,flag);
            %
            % INPUT
            % ------------
            %
            % A    : the input nonnegative m x n matrix A
            % k    : the rank of the computed factors W,H
            % flag : indicates the variant of the NNDSVD Algorithm
            %        flag = 0 --> NNDSVD
            %        flag = 1 --> NNDSVDa
            %        flag = 2 --> NNDSVDar
            %        flag = 3 --> NNDSVD using random SVD calculation
            %
            % OUTPUT
            % -------------
            %
            % W   : nonnegative m x k matrix
            % H   : nonnegative k x n matrix
            %
            %
            % References:
            %
            % [1] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
            %     start for nonnegative matrix factorization, Pattern Recognition,
            %     Elsevier
            %
            % This code is kindly provided by the authors for research porpuses.
            % - Efstratios Gallopoulos (stratis@ceid.upatras.gr)
            % - Christos Boutsidis (boutsc@cs.rpi.edu)
            %
            % For any problems or questions please send an email to boutsc@cs.rpi.edu
            %--------------------------------------------------------------------------

            import mladni.NMFManopt.pos
            import mladni.NMFManopt.neg

            %----------------------check the input matrix------------------------------
            if numel(find(A<0)) > 0
                error('The input matrix contains negative elements !')
            end
            %--------------------------------------------------------------------------

            %size of the input matrix
            [m,n] = size(A);

            %the matrices of the factorization
            W = zeros(m,k);
            H = zeros(k,n);

            % 1st SVD --> partial SVD rank-k to the input matrix A.
            if ( flag == 3)
                % use random svd for efficient computation
                l = max(3*k,20) ;
                [U,S,V]=randpca(A,k,true,8,l);
            else
                % use standard matlab svn implementation
                [U,S,V] = svds(A,k);
            end

            %choose the first singular triplet to be nonnegative
            W(:,1)     =  sqrt(S(1,1)) * abs(U(:,1) );
            H(1,:)     =  sqrt(S(1,1)) * abs(V(:,1)');

            % 2nd SVD for the other factors (see table 1 in our paper)
            for i=2:k
                uu = U(:,i); vv = V(:,i);
                uup = pos(uu); uun = neg(uu) ;
                vvp = pos(vv); vvn = neg(vv);
                n_uup = norm(uup);
                n_vvp = norm(vvp) ;
                n_uun = norm(uun) ;
                n_vvn = norm(vvn) ;
                termp = n_uup*n_vvp; termn = n_uun*n_vvn;
                if (termp >= termn)
                    W(:,i) = sqrt(S(i,i)*termp)*uup/n_uup;
                    H(i,:) = sqrt(S(i,i)*termp)*vvp'/n_vvp;
                else
                    W(:,i) = sqrt(S(i,i)*termn)*uun/n_uun;
                    H(i,:) = sqrt(S(i,i)*termn)*vvn'/n_vvn;
                end
            end
            %------------------------------------------------------------

            %actually these numbers are zeros
            W(W<0.0000000001)=0.1;
            H(H<0.0000000001)=0.1;

            if(exist('flag','var'))
                % NNDSVDa: fill in the zero elements with the average
                if flag==1
                    ind1      =  W==0 ;
                    ind2      =  H==0 ;
                    average   =  mean(A(:)) ;
                    W( ind1 ) =  average    ;
                    H( ind2 ) =  average    ;

                    % NNDSVDar: fill in the zero elements with random values in the space [0:average/100]
                elseif flag==2
                    ind1      =  find(W==0) ;
                    ind2      =  find(H==0) ;
                    n1        =  numel(ind1);
                    n2        =  numel(ind2);

                    average   =  mean(A(:))       ;
                    W( ind1 ) =  (average*rand(n1,1)./100)  ;
                    H( ind2 ) =  (average*rand(n2,1)./100)  ;
                end
            end
        end

        function [Ap] = pos(A)
            %% This function sets to zero the negative elements of a matrix
            Ap = (A>=0).*A;
        end

        function [Am] = neg(A)
            %% This functions sets to zero the positive elements of a matrix and takes
            %  the absolute value of the negative elements
            Am = (A<0).*(-A);
        end

        function M = positivegrassmannfactory(n, p, k, gpuflag)
            % Returns a manifold struct to optimize over the space of vector subspaces,
            % the open set of the positive orthants.
            %
            % function M = positivegrassmannfactory(n, p)
            % function M = positivegrassmannfactory(n, p, k)
            % function M = positivegrassmannfactory(n, p, k, gpuflag)
            %
            % Grassmann manifold: each point on this manifold is a collection of k
            % vector subspaces of dimension p embedded in R^n.
            %
            % The metric is obtained by making the Grassmannian a Riemannian quotient
            % manifold of the Stiefel manifold, i.e., the manifold of orthonormal
            % matrices, itself endowed with a metric by making it a Riemannian
            % submanifold of the Euclidean space, endowed with the usual inner product.
            % In short: it is the usual metric used in most cases.
            %
            % This structure deals with matrices X of size n x p x k (or n x p if
            % k = 1, which is the default) such that each n x p matrix is orthonormal,
            % i.e., X'*X = eye(p) if k = 1, or X(:, :, i)' * X(:, :, i) = eye(p) for
            % i = 1 : k if k > 1. Each n x p matrix is a numerical representation of
            % the vector subspace its columns span.
            %
            % The retraction is based on a polar factorization and is second order.
            %
            % Set gpuflag = true to have points, tangent vectors and ambient vectors
            % stored on the GPU. If so, computations can be done on the GPU directly.
            %
            % By default, k = 1 and gpuflag = false.
            %
            % See also: stiefelfactory grassmanncomplexfactory grassmanngeneralizedfactory

            % This file is part of Manopt: www.manopt.org.
            % Original author: Nicolas Boumal, Dec. 30, 2012.
            % Contributors:
            % Change log:
            %   March 22, 2013 (NB):
            %       Implemented geodesic distance.
            %
            %   April 17, 2013 (NB):
            %       Retraction changed to the polar decomposition, so that the vector
            %       transport is now correct, in the sense that it is compatible with
            %       the retraction, i.e., transporting a tangent vector G from U to V
            %       where V = Retr(U, H) will give Z, and transporting GQ from UQ to VQ
            %       will give ZQ: there is no dependence on the representation, which
            %       is as it should be. Notice that the polar factorization requires an
            %       SVD whereas the qfactor retraction requires a QR decomposition,
            %       which is cheaper. Hence, if the retraction happens to be a
            %       bottleneck in your application and you are not using vector
            %       transports, you may want to replace the retraction with a qfactor.
            %
            %   July  4, 2013 (NB):
            %       Added support for the logarithmic map 'log'.
            %
            %   July  5, 2013 (NB):
            %       Added support for ehess2rhess.
            %
            %   June 24, 2014 (NB):
            %       Small bug fix in the retraction, and added final
            %       re-orthonormalization at the end of the exponential map. This
            %       follows discussions on the forum where it appeared there is a
            %       significant loss in orthonormality without that extra step. Also
            %       changed the randvec function so that it now returns a globally
            %       normalized vector, not a vector where each component is normalized
            %       (this only matters if k>1).
            %
            %   July 8, 2018 (NB):
            %       Inverse retraction implemented.
            %
            %   Aug. 3, 2018 (NB):
            %       Added GPU support: just set gpuflag = true.
            %
            %   Apr. 19, 2019 (NB):
            %       ehess2rhess: to ensure horizontality, it makes sense to project
            %       last, same as in stiefelfactory.
            %
            %   May 3, 2019 (NB):
            %       Added explanation about vector transport relation to retraction.
            %
            %   Nov. 13, 2019 (NB):
            %       Added pairmean function.
            %
            %   Jan. 8, 2021 (NB)
            %       Added tangent2ambient/tangent2ambient_is_identity pair.
            %       Here, 'ambient' refers to the total space.
            %
            %   May 14, 2024 (JJ Lee)
            %       Comparing to positivefactory and experimenting, M.{inner,norm,proj,retraction,transp}  
            %       appear to be the objects requiring adjustments to maintain positivity.   
            %       I've tried NB's Riemannian metric from Intro. Optimization Manifolds, sec 11.6,
            %       which risks division by zero, but this fails to solve because of 
            %       negative curvatures.

            assert(n >= p, ...
                ['The dimension n of the ambient space must be larger ' ...
                'than the dimension p of the subspaces.']);

            if ~exist('k', 'var') || isempty(k)
                k = 1;
            end
            if ~exist('gpuflag', 'var') || isempty(gpuflag)
                gpuflag = false;
            end

            % If gpuflag is active, new arrays (e.g., via rand, randn, zeros, ones)
            % are created directly on the GPU; otherwise, they are created in the
            % usual way (in double precision).
            if gpuflag
                array_type = 'gpuArray';
            else
                array_type = 'double';
            end

            if k == 1
                M.name = @() sprintf('positive Grassmann manifold Gr(%d, %d)', n, p);
            elseif k > 1
                M.name = @() sprintf('positive Multi Grassmann manifold Gr(%d, %d)^%d', ...
                    n, p, k);
            else
                error('k must be an integer no less than 1.');
            end

            M.dim = @() k*p*(n-p);

            %%% M.inner = @(x, d1, d2) d1(:).'*d2(:);
            M.inner = @pinner;
            function pin = pinner(x, d1, d2)
                D1 = exp(log(d1) - log(x));
                D2 = exp(log(d2) - log(x));
                pin = D1(:).'*D2(:);
            end

            %%% M.norm = @(x, d) norm(d(:));
            M.norm = @(x, d) sqrt(pinner(x, d, d)); 

            M.dist = @distance;
            function d = distance(x, y)
                square_d = 0;
                XtY = multiprod(multitransp(x), y);
                for kk = 1 : k
                    cos_princ_angle = svd(XtY(:, :, kk));
                    % For x and y closer than ~sqrt(eps), this function is
                    % inaccurate, and typically returns values close to ~sqrt(eps).
                    square_d = square_d + sum(real(acos(cos_princ_angle)).^2);
                end
                %assert(~any(square_d < 0, "all"))  % satisfied
                d = sqrt(square_d);
            end

            M.typicaldist = @() sqrt(p*k);

            % Orthogonal projection of an ambient vector U to the horizontal space
            % at X.
            M.proj = @projection;
            function Up = projection(X, U)

                XtU = multiprod(multitransp(X), U);
                Up = U - multiprod(X, XtU);
                %assert(~any(Up < 0, "all"))  % not satisfied
            end

            M.tangent = M.proj;

            M.tangent2ambient_is_identity = true;
            M.tangent2ambient = @(X, U) U;

            M.egrad2rgrad = M.proj;

            M.ehess2rhess = @ehess2rhess;
            function rhess = ehess2rhess(X, egrad, ehess, H)
                XtG = multiprod(multitransp(X), egrad);
                HXtG = multiprod(H, XtG);
                rhess = projection(X, ehess - HXtG);
            end

            M.retr = @retraction;
            function Y = retraction(X, U, t)
                if nargin < 3
                    Y = X + U;
                else
                    Y = X + t*U;
                end
                for kk = 1 : k

                    % Compute the polar factorization of Y = X+tU
                    [u, s, v] = svd(Y(:, :, kk), 'econ'); %#ok
                    Y(:, :, kk) = u*v';

                    % Another way to compute this retraction uses QR instead of SVD.
                    % As compared with the Stiefel factory, we do not need to
                    % worry about flipping signs of columns here, since only
                    % the column space is important, not the actual columns.
                    % We prefer the polar factor to the Q-factor computation for
                    % reasons explained below: see M.transp.
                    %
                    % [Q, unused] = qr(Y(:, :, kk), 0); %#ok
                    % Y(:, :, kk) = Q;

                end
                % assert(~any(Y < 0, "all"))  % not satisfied
            end

            % This inverse retraction is valid for both the QR retraction and the
            % polar retraction.
            M.invretr = @invretr;
            function U = invretr(X, Y)
                XtY = multiprod(multitransp(X), Y);
                U = zeros(n, p, k, array_type);
                for kk = 1 : k
                    U(:, :, kk) = Y(:, :, kk) / XtY(:, :, kk);
                end
                U = U - X;
                % assert(~any(Up < 0, "all"))  % satisfied
            end

            % See Eq. (2.65) in Edelman, Arias and Smith 1998.
            M.exp = @exponential;
            function Y = exponential(X, U, t)
                if nargin == 3
                    tU = t*U;
                else
                    tU = U;
                end
                Y = zeros(size(X), array_type);
                for kk = 1 : k
                    [u, s, v] = svd(tU(:, :, kk), 0);
                    cos_s = diag(cos(diag(s)));
                    sin_s = diag(sin(diag(s)));
                    Y(:, :, kk) = X(:, :, kk)*v*cos_s*v' + u*sin_s*v';
                    % From numerical experiments, it seems necessary to
                    % re-orthonormalize. This is overall quite expensive.
                    [q, unused] = qr(Y(:, :, kk), 0); %#ok
                    Y(:, :, kk) = q;
                end
            end

            % Test code for the logarithm:
            % Gr = grassmannfactory(5, 2, 3);
            % x = Gr.rand()
            % y = Gr.rand()
            % u = Gr.log(x, y)
            % Gr.dist(x, y) % These two numbers should
            % Gr.norm(x, u) % be the same.
            % z = Gr.exp(x, u) % z needs not be the same matrix as y, but it should
            % v = Gr.log(x, z) % be the same point as y on Grassmann: dist almost 0.
            M.log = @logarithm;
            function U = logarithm(X, Y)
                U = zeros(n, p, k, array_type);
                for kk = 1 : k
                    x = X(:, :, kk);
                    y = Y(:, :, kk);
                    ytx = y.'*x;
                    At = y.'-ytx*x.';
                    Bt = ytx\At;
                    [u, s, v] = svd(Bt.', 'econ');

                    u = u(:, 1:p);
                    s = diag(s);
                    s = s(1:p);
                    v = v(:, 1:p);

                    U(:, :, kk) = u*diag(atan(s))*v.';
                end
            end

            M.hash = @(X) ['z' hashmd5(X(:))];

            M.rand = @prandom;
            function X = random()
                X = randn(n, p, k, array_type);
                for kk = 1 : k
                    [Q, unused] = qr(X(:, :, kk), 0); %#ok
                    X(:, :, kk) = Q;
                end
            end
            function X = prandom()
                X = exp(randn(n, p, k, array_type));
                for kk = 1 : k
                    [Q, unused] = qr(X(:, :, kk), 0); %#ok
                    X(:, :, kk) = Q;
                end
                X(X < eps) = eps;
            end

            M.randvec = @prandomvec;
            function U = randomvec(X)
                U = projection(X, randn(n, p, k, array_type));
                U = U / norm(U(:));
            end
            function U = prandomvec(X)
                U = projection(X, randn(n, p, k, array_type));
                U = U / M.norm(X, U);
                U(U < eps) = eps;
            end

            M.lincomb = @matrixlincomb;

            M.zerovec = @(x) zeros(n, p, k, array_type);

            % This transport is compatible with the polar retraction, in the
            % following sense:
            %
            % n = 7; p = 3;
            % Gr = grassmannfactory(n, p);
            % X = Gr.rand();
            % U = Gr.randvec(X);
            % V = Gr.randvec(X);
            % [Q, ~] = qr(randn(p));
            % Gr.transp(X*Q, Gr.retr(X*Q, V*Q), U*Q) % these two
            % Gr.transp(X, Gr.retr(X, V), U)*Q       % are equal (up to eps)
            %
            % That is, if we transport U, the horizontal lift of some tangent
            % vector at X, to Y, and Y = Retr_X(V) with V the horizontal lift of
            % some tangent vector at X, we get the horizontal lift of some tangent
            % vector at Y. If we displace X, U, V to XQ, UQ, VQ for some arbitrary
            % orthogonal matrix Q, we get a horizontal lift of some vector at YQ.
            % Importantly, these two vectors are the lifts of the same tangent
            % vector, only lifted at Y and YQ.
            %
            % However, this vector transport is /not/ fully invariant, in the sense
            % that transporting U from X to some arbitrary Y may well yield the
            % lift of a different vector when compared to transporting U from X
            % to YQ, where Q is an arbitrary orthogonal matrix, even though YQ is
            % equivalent to Y. Specifically:
            %
            % Y = Gr.rand();
            % Gr.transp(X, Y*Q, U) - Gr.transp(X, Y, U)*Q   % this is not zero.
            %
            % However, the following vectors are equal:
            %
            % Gr.transp(X, Y*Q, U) - Gr.transp(X, Y, U)     % this *is* zero.
            %
            % For this to be a proper vector transport from [X] to [Y] in general,
            % assuming X'Y is invertible, one should multiply the output of this
            % function on the right with the polar factor of X'*Y, that is,
            % multiply by u*v' where [u, s, v] = svd(X'*Y), for each slice.
            M.transp = @(X, Y, U) projection(Y, U);

            % The mean of two points is here defined as the midpoint of a
            % minimizing geodesic connecting the two points. If the log of (X1, X2)
            % is not uniquely defined, then the returned object may not be
            % meaningful; in other words: this works best if (X1, X2) are close.
            M.pairmean = @pairmean;
            function Y = pairmean(X1, X2)
                Y = M.exp(X1, .5*M.log(X1, X2));
            end

            M.vec = @(x, u_mat) u_mat(:);
            M.mat = @(x, u_vec) reshape(u_vec, [n, p, k]);
            M.vecmatareisometries = @() true;


            % Automatically convert a number of tools to support GPU.
            if gpuflag
                M = factorygpuhelper(M);
            end
        end
            
        function [X_, r_] = rand_X(opts)
            %% generate nonnegative random data for illustration purposes.

            arguments
                opts.m {mustBeScalarOrEmpty} = 42
                opts.n {mustBeScalarOrEmpty} = 60
                opts.r {mustBeScalarOrEmpty} = 5
            end

            X_ = randn(opts.m, opts.n);
            r_ = opts.r;
        end        
    end


    %% PRIVATE

    methods (Access = private)
        function [W, info_] = nmf(this, X, r, do_AD, epsilon)
            % Attempts nonnegative matrix factorization of X truncated to rank r.
            %
            % function [W, info] = nmf(X, r, do_AD, epsilon)
            %
            % Input: A real matrix X of size mxn and an integer r <= min(m, n).
            %        Logical do_AD implements manoptAD().  
            %        Real epsilon provides smoothing:  try 0 or 1.
            % Output: An orthonormal matrix W of size mxr, such that W*W'*X is the best rank r
            %         approximation of X according to the Frobenius norm. All real.
            %
            % The decomposition is obtained by minimizing
            %   f(W) = norm(W*W' - X, 'fro')^2
            % with W orthonormal. Notice that f(W*Q) = f(W) for all
            % Q orthogonal pxp matrices. Hence, only the column space of W
            % matters and we may perform the optimization over a Grassmannian manifold. 
            % We are going for a best low-rank approximation of X.
            %
            % The inner workings of the Grassmann manifold use the built-in svd
            % function of Matlab but only for matrices of size mxp and nxp to
            % re-orthonormalize them.
            %
            % Notice that we are actually chasing a best fixed-rank approximation of a
            % matrix, which is best obtained by working directly over a manifold of
            % fixed-rank matrices. This is simply an example script to demonstrate some
            % functionalities of the toolbox.
            %
            % The code can be modified to accept a function handle for X(x) = X*x
            % instead of a matrix X, which is often useful. This would further require
            % a function handle Xt for the transpose of X, such that Xt(x) = X.'*x.

            % The original file is part of Manopt and is copyrighted. See the license file.
            %
            % Main author: Nicolas Boumal, July 5, 2013
            % Contributors: John J. Lee, May 14, 2024
            %
            % Change log:
            %
            %   Aug. 23, 2021 (XJ):
            %       Added AD to compute the egrad and the ehess
            %   May 14, 2024 (JJL):
            %       Adjusted for trials of NMF

            assert(~any(X < 0, "all"))
            do_force_nonneg_ = this.do_force_nonneg;
            do_yang_oja_ = this.do_yang_oja;  % avoid exposing this to inner functions

            % Retrieve the size of the problem and make sure the requested
            % approximation rank is at most the maximum possible rank.
            [m_, n_] = size(X);
            assert(r <= min(m_, n_), 'r must be smaller than the smallest dimension of X.');

            % Define a problem structure, specifying the manifold M, the cost
            % function and its derivatives. Here, to demonstrate the rapid
            % prototyping capabilities of Manopt, we directly define the Euclidean
            % gradient and the Euclidean Hessian egrad and ehess instead of the
            % Riemannian gradient and Hessian grad and hess. Manopt will take care
            % of the conversion. This automatic conversion is usually not
            % computationally optimal though, because much of the computations
            % involved in obtaining the gradient could be reused to obtain the
            % Hessian. After the prototyping stage, when efficiency becomes
            % important, it makes sense to define grad and hess rather than egrad
            % an ehess, and to use the caching system (the store structure).

            % All of the code will work just as well if we ignore the invariance
            % property of the cost function indicated above and thus place W
            % on the Stiefel manifold (orthonormal matrices) instead of the
            % Grassmann manifold. Working on Stiefel is expected to be slower
            % though, partly because de search space is higher dimensional and
            % partly because the optimizers are not isolated.
            % problem_.M = stiefelfactory(m_, r);
            problem_.M = this.manifold_factory(m_, r);

            if do_AD

                % An alternative way to compute the egrad and the ehess is to use
                % automatic differentiation provided in the deep learning toolbox
                % (slower). Notice that the function norm is not supported for AD so
                % far. Replace norm(...,'fro') with the backup function cnormsqfro
                % described in manoptADhelp
                problem_.cost = @cost_AD;
                
                % call manoptAD to prepare AD for the problem structure
                problem_ = manoptAD(problem_, 'nohess');
            else
                problem_.cost  = @cost;
                problem_.egrad = @egrad;
            end

            if this.do_plot
                this.problem = problem_;
                plot_checks(this)
            end

            % Issue a call to a solver. A random initial guess will be chosen and
            % default options are selected. Here, we specify a maximum trust
            % region radius (which in turn induces an initial trust region radius).
            % Note that this is not required: default values are used if we omit
            % this. The diameter of the manifold scales like sqrt(2*r), hence the
            % form of our (empirical) choice.
            options.Delta_bar = 1e5*sqrt(2*r);  % 4 -> 1e5
            options.tolgradnorm = 1e-7;
        	options.verbosity = 3; % Change this number for more or less output

            % Starting guess for W uses NNDSVD using random SVD calculation
            % https://github.com/asotiras/brainparts/blob/master/NNDSVD.m
            W = this.W0;
            if epsilon > 0
                n_iter = this.n_iterations;
            else
                n_iter = 1;
            end
            for iter = 1:n_iter
                switch char(this.solver_name)
                    case 'trustregions'
                        options.maxiter = 1e3;
                        [W, ~, info_] = trustregions(problem_, W, options);
                    case 'arc'
                        options.maxiter = 1e3;
                        [W, ~, info_] = arc(problem_, W, options);  
                        % The subproblem solver failed to make progress even on the model; this is likely due to numerical errors.
                    case 'conjugategradient'
                        options.maxiter = 1e3;
                        [W, ~, info_] = conjugategradient(problem_, W, options);
                    case 'pso'
                        [W, ~, info_] = pso(problem_, W, options);
                    case 'stochasticgradient'
                        %% see also PCA_stochastic
                        options.checkperiod = min(1e3, floor(this.maxiter/10));
                        options.statsfun = statsfunhelper('metric', @cost);
                        options.maxiter = this.maxiter;
                        options.batchsize = this.batchsize;
                        options.stepsize_type = 'decay';
                        options.stepsize_init = 1e2;
                        options.stepsize_lambda = this.stepsize_lambda;
                        problem_.ncostterms = n_;
                        problem_.cost = [];
                        problem_.partialegrad = @partialegrad;
                        [W, info_] = stochasticgradient(problem_, W, options);
                    otherwise
                        error("mladni:ValueError", stackstr())
                end
                epsilon = epsilon * this.reduction;
            end

            %% collect results for this

            this.problem = problem_;
            this.info = info_;
            this.result = W;
            epsilon = 0;
            this.cost_result = cost(W);

            %% INNER FUNCTIONS

            % The functions below make many redundant computations. This
            % performance hit can be alleviated by using the caching system.

            % Cost function
            function f = cost(W)
                % W(W<eps)=eps;  % breaks gradients
                WtX = W.'*X;  % r x n
                vecs = W*WtX - X;  % m x n
                sqnrms = sum(vecs.^2, 1);
                vals = sqrt(sqnrms + epsilon^2) - epsilon;
                f = mean(vals);
            end

            % Cost function with auto-differentiation
            function f = cost_AD(W)
                % W(W<eps)=eps;  % breaks gradients
                WtX = W.'*X;  % r x n
                f = cnormsqfro(W*WtX - X) + epsilon^2;  % m x n
            end

            % Euclidean gradient of the cost function
            function G = egrad(W)

                %W(W<eps)=eps;

                if do_yang_oja_
                    % https://github.com/asotiras/brainparts/blob/master/opnmf.m
                    XtW = X.'*W;  % n x r
                    numer = X*XtW;  % m x r
                    Wtnumer = W.'*numer;  % r x r
                    denom = W*Wtnumer;  % m x r
                    update = numer./denom;  % m x r
                    % multiplicative update rule
                    W = W.*update;
                    % As the iterations were progressing, computational time per iteration was increasing due to operations involving really small values
                    W(W<eps)=eps;
                    W = W./norm(W, 'fro');
                end

        		% Note that the computation of vecs and sqnrms is redundant
        		% with their computation in the cost function. To speed
        		% up the code, it would be wise to use the caching capabilities
        		% of Manopt (the store structure). See online documentation.
        		% It is not done here to keep the code a bit simpler.

                % W = abs(W);  % produces negative curvature

                WtX = W.'*X;  % r x n
                vecs = W*WtX - X;  % m x n; sign reversal produces negative curvature
                sqnrms = sum(vecs.^2, 1); 

                % This explicit loop is a bit slow: the code below is equivalent
                % and faster to compute the gradient.
                % G = zeros(m, r);
                % for i=1:n
                %     G = G + (1/sqrt(sqnrms(i) + epsilon^2)) * vecs(:,i) * WtX(:,i)';
                % end
                % G = G/n;
                G = mean( ...
                    multiscale(1./sqrt(sqnrms + epsilon^2), ...
                    multiprod(reshape(vecs, [m_, 1, n_]), ...
                    multitransp(reshape(WtX, [r, 1, n_])))), 3);
            end

            % Euclidean gradient of the cost function for sample of partial data
            function G = partialegrad(W, sample)
                % W is an orthonormal matrix of size mxr.
                % Sample is a vector of indices between 1 and n: a subset.
                % Extract a subset of the dataset X.
                
                % Forcing nonnegativity is tolerable for SGD, 
                % but other gradient methods, that also use consistently defined costs, fail gradient checks.
                if do_force_nonneg_
                    W(W<eps)=eps;
                end

                X_ = X(:, sample);
                nsample = size(X_, 2);

                if do_yang_oja_
                    XtW = X.'*W;  % n x r
                    numer = X*XtW;  % m x r
                    Wtnumer = W.'*numer;  % r x r
                    denom = W*Wtnumer;  % m x r
                    update = numer./denom;  % m x r
                    % multiplicative update rule
                    W = W.*update;
                    % As the iterations were progressing, computational time per iteration was increasing due to operations involving really small values
                    W(W<eps)=eps;
                    W = W./norm(W, 'fro');
                end

                WtX = W.'*X_;  % r x n
                vecs = W*WtX - X_;  % m x n; sign reversal produces negative curvature
                sqnrms = sum(vecs.^2, 1); 

                G = mean( ...
                    multiscale(1./sqrt(sqnrms + epsilon^2), ...
                    multiprod(reshape(vecs, [m_, 1, nsample]), ...
                    multitransp(reshape(WtX, [r, 1, nsample])))), 3);
            end
        end

        function [U, S, V, info_] = truncated_svd(this, A, p, do_AD)
            % Returns an SVD decomposition of A truncated to rank p.
            %
            % function [U, S, V, info] = truncated_svd(A, p, do_AD)
            %
            % Input: A real matrix A of size mxn and an integer p <= min(m, n).
            %        Logical do_AD implements manoptAD().  
            % Output: An orthonormal matrix U of size mxp, an orthonormal matrix Y of
            %         size nxp and a diagonal matrix S of size pxp with nonnegative and
            %         decreasing diagonal entries such that USV.' is the best rank p
            %         approximation of A according to the Frobenius norm. All real.
            %         This function produces an output akin to svds.
            %
            % The decomposition is obtained by maximizing
            %   f(U, V) = .5*norm(U'*A*V, 'fro')^2
            % where U, V are orthonormal. Notice that f(U*Q, V*R) = f(U, V) for all
            % Q, R orthogonal pxp matrices. Hence, only the column spaces of U and V
            % matter and we may perform the optimization over a product of two
            % Grassmannian manifolds.
            %
            % It is easy to show that maximizing f is equivalent to minimizing g with
            %   g(U, V) = min_S norm(U*S*V' - A, 'fro')^2,
            % which confirms that we are going for a best low-rank approximation of A.
            %
            % The inner workings of the Grassmann manifold use the built-in svd
            % function of Matlab but only for matrices of size mxp and nxp to
            % re-orthonormalize them.
            %
            % Notice that we are actually chasing a best fixed-rank approximation of a
            % matrix, which is best obtained by working directly over a manifold of
            % fixed-rank matrices. This is simply an example script to demonstrate some
            % functionalities of the toolbox.
            %
            % The code can be modified to accept a function handle for A(x) = A*x
            % instead of a matrix A, which is often useful. This would further require
            % a function handle At for the transpose of A, such that At(x) = A.'*x.

            % This file is part of Manopt and is copyrighted. See the license file.
            %
            % Main author: Nicolas Boumal, July 5, 2013
            % Contributors:
            %
            % Change log:
            %
            %   Aug. 23, 2021 (XJ):
            %       Added AD to compute the egrad and the ehess

            % Retrieve the size of the problem and make sure the requested
            % approximation rank is at most the maximum possible rank.
            [m_, n_] = size(A);
            assert(p <= min(m_, n_), 'p must be smaller than the smallest dimension of A.');

            % Define the cost and its derivatives on the Grassmann manifold
            tuple.U = this.manifold_factory(m_, p);
            tuple.V = this.manifold_factory(n_, p);

            % All of the code will work just as well if we ignore the invariance
            % property of the cost function indicated above and thus place U and V
            % on the Stiefel manifold (orthonormal matrices) instead of the
            % Grassmann manifold. Working on Stiefel is expected to be slower
            % though, partly because de search space is higher dimensional and
            % partly because the optimizers are not isolated.
            % tuple.U = stiefelfactory(m, p);
            % tuple.V = stiefelfactory(n, p);
            M = productmanifold(tuple);

            % Define a problem structure, specifying the manifold M, the cost
            % function and its derivatives. Here, to demonstrate the rapid
            % prototyping capabilities of Manopt, we directly define the Euclidean
            % gradient and the Euclidean Hessian egrad and ehess instead of the
            % Riemannian gradient and Hessian grad and hess. Manopt will take care
            % of the conversion. This automatic conversion is usually not
            % computationally optimal though, because much of the computations
            % involved in obtaining the gradient could be reused to obtain the
            % Hessian. After the prototyping stage, when efficiency becomes
            % important, it makes sense to define grad and hess rather than egrad
            % an ehess, and to use the caching system (the store structure).
            problem_.M = M;

            if do_AD

                % An alternative way to compute the egrad and the ehess is to use
                % automatic differentiation provided in the deep learning toolbox
                % (slower). Notice that the function norm is not supported for AD so
                % far. Replace norm(...,'fro') with the backup function cnormsqfro
                % described in manoptADhelp
                problem_.cost = @cost_AD;
                
                % call manoptAD to prepare AD for the problem structure
                problem_ = manoptAD(problem_);
            else
                problem_.cost  = @cost;
                problem_.egrad = @egrad;
                problem_.ehess = @ehess;
            end

            if this.do_plot
                this.problem = problem_;
                plot_checks(this)
            end

            % Issue a call to a solver. A random initial guess will be chosen and
            % default options are selected. Here, we specify a maximum trust
            % region radius (which in turn induces an initial trust region radius).
            % Note that this is not required: default values are used if we omit
            % this. The diameter of the manifold scales like sqrt(2*p), hence the
            % form of our (empirical) choice.
            options.Delta_bar = 4*sqrt(2*p);
            % options.tolgradnorm = 1e-7;
        	options.verbosity = 3; % Change this number for more or less output
            % options.maxiter = 1e4;
            X_cost_ = [];
            % [X_, X_cost_, info_] = trustregions(problem_, [], options);
            [X_, X_cost_, info_] = arc(problem_, [], options);
            % [X_, X_cost_, info_] = conjugategradient(problem_, [], options);
            % [X_, info_] = stochasticgradient(problem_, [], options);
            % [X_, X_cost_, info_] = pso(problem_, [], options);
            U = X_.U;
            V = X_.V;

            % Finish the job by rotating U and V such that the middle matrix S can
            % be diagonal with nonnegative, decreasing entries. This requires a
            % small svd of size pxp.
            Spp = U'*A*V;
            [Upp, Spp, Vpp] = svd(Spp);
            U = U*Upp;
            S = Spp;
            V = V*Vpp;

            %% collect results for this

            this.problem = problem_;
            this.info = info_;
            this.result = X_;
            this.cost_result = X_cost_;

            %% INNER FUNCTIONS

            % The functions below make many redundant computations. This
            % performance hit can be alleviated by using the caching system.

            % Cost function
            function f = cost(X__)
                f = -.5*norm(X__.U'*A*X__.V, 'fro')^2;
            end

            % Cost function with auto-differentiation
            function f = cost_AD(X__)
                f = -.5*cnormsqfro(X__.U'*A*X__.V);
            end

            % Euclidean gradient of the cost function
            function g = egrad(X__)
                U__ = X__.U;
                V__ = X__.V;
                AV = A*V__;
                AtU = A'*U__;
                g.U = -AV*(AV'*U__);
                g.V = -AtU*(AtU'*V__);
            end

            % Euclidean Hessian of the cost function
            function h = ehess(X__, H__)
                U__ = X__.U;
                V__ = X__.V;
                Udot = H__.U;
                Vdot = H__.V;
                AV = A*V__;
                AtU = A'*U__;
                AVdot = A*Vdot;
                AtUdot = A'*Udot;
                h.U = -(AVdot*AV'*U__ + AV*AVdot'*U__ + AV*AV'*Udot);
                h.V = -(AtUdot*AtU'*V__ + AtU*AtUdot'*V__ + AtU*AtU'*Vdot);
            end            

            % In this helper function, given a point 'X' on the manifold we check
            % whether the caching structure 'store' has been populated with
            % quantities that are useful to compute at X or not. If they were not,
            % then we compute and store them now.
            function store = prepare(X__, store)
                if ~isfield(store, 'ready') || ~store.ready
                    store.Ut_A_V = X__.U'*A*X__.V;
                    store.ready = true;
                end
            end
        end        
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
