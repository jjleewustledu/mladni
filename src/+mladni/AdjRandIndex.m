classdef AdjRandIndex
    %% line1
    %  line2
    %  
    %  Created 09-Jan-2024 21:42:39 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        dataset_home
    end

    methods
        function this = AdjRandIndex(dataset_home)
            arguments
                dataset_home {mustBeFolder} = pwd
            end
            this.dataset_home = dataset_home;
        end

        function clusterdata(this)
        end
        function pdist(this)
        end
        function linkage(this)
        end
        function cluster(this)
        end
        
    end

    methods (Static)
        function r = clustering_fast(u,v)
            % clustering quality measures assumptions : 
            % 0 corresponds to background
            % we do not care about the background
            % cluster labels are assumed to be enumerated from 1 to max number of
            % clusters
            %
            % this function should not be used when comparing binary segmentations
            
            m=max(max(u),max(v));
            
            if(m == 1)
                error('AdjRandIndex:argChk','clustering_fast() should not be used for comparing binary segmentations');
            end
            
            va=zeros(1,m);
            vb=zeros(1,m);
            mat=zeros(m);
            
            for i=1:m
                va(i) = sum(u==i) ;
                vb(i) = sum(v==i) ;
                hu = (u==i) ;
                for j=1:m
                    hv = (v==j) ;
                    mat(i,j) = sum(hu.*hv);
                end
            end
            
            ra=sum(sum(mat.*(mat-ones(m))))/2.0;
            
            rb=va*(va-ones(1,m))'/2.0;
            
            rc=vb*(vb-ones(1,m))'/2.0;
            
            rn=length(u)*(length(u)-1)/2.0;
            
            r=(ra-rb*rc/rn)/( 0.5*rb+0.5*rc-rb*rc/rn );
        end  
        
        function [ARI,overlap] = evaluate_pair(B1, B2, opts)
            %% evaluates the adjusted Rand Index (ARI) between NMF results B1 and B2, equally shaped.
            %
            %  ARI is evaluated by deriving hard clusters from the estimated components.
            %
            %  Args:
            %    B1 {mustBeNumeric}  % Nvoxels x Ncomponents; Ncomponents may be 1
            %    B2 {mustBeNumeric}  % same shape as B1
            %    opts.do_hungarian logical = true
            %    opts.do_ari logical = true            
            % 
            %  Returns:
            %    ARI : scalar
            %    overlap : Ncomponents x 1

            arguments
                B1 {mustBeNumeric}  % Nvoxels x Ncomponents
                B2 {mustBeNumeric}  % Nvoxels x Ncomponents
                opts.do_hungarian logical = true
                opts.do_ari logical = true
            end
            
            % normalize to unit norm
            wlen1 = sqrt(sum((B1).^2)) ;  % 1 x Ncomponents
            wlen2 = sqrt(sum((B2).^2)) ;  % same shape as wlen1
            
            if any(wlen1==0)
                wlen1(wlen1==0) = 1;
            end
            W1 = bsxfun(@times,B1,1./wlen1) ;  % normalize
           
            if any(wlen2==0)
                wlen2(wlen2==0) = 1;
            end
            W2 = bsxfun(@times,B2,1./wlen2) ;  % normalize
            
            if opts.do_hungarian
                % calculate inner products
                inner_product = W1'*W2 ;  % Ncomponents x Ncomponents

                % take a distance
                dist = 2*(1 - inner_product) ;

                % find correspondences
                [Matching,~] = mladni.AdjRandIndex.Hungarian(dist);
                [~,idx_hug1]=max(Matching,[],2);

                % overlap - hungarian
                overlap = zeros(length(wlen1),1) ;
                for b=1:length(wlen1)
                    overlap(b) = inner_product(b,idx_hug1(b));
                end

                % overlap with best
                %overlap_best = max(inner_product,[],2) ;
            else
                overlap = nan;
            end
            
            % also evaluate overlap based on adjusted Rand Index    
            rowLen1 = sum(W1,2) ;
            rowLen2 = sum(W2,2) ;
            
            if any(rowLen1==0)
                rowLen1(rowLen1==0) = 1 ;
            end
            if any(rowLen2==0)
                rowLen2(rowLen2==0) = 1 ;
            end
            WW1 = bsxfun(@times,(W1'),1./(rowLen1')); WW1=WW1';
            WW2 = bsxfun(@times,(W2'),1./(rowLen2')); WW2=WW2';
            
            [~,clustering1] = max(WW1,[],2);
            [~,clustering2] = max(WW2,[],2);
            ARI = mladni.AdjRandIndex.clustering_fast(clustering1,clustering2);   
        end

        function S = inspect_bimodal_ARI(opts)
            %% provides more information for revealing why ARI reproducibility appear bimodal after bootstrapping

            arguments
                opts.span {mustBeNumeric} = 24
                opts.anticlust_home {mustBeFolder} = ...
                    fullfile(getenv("SINGULARITY_HOME"), "ADNI", "NMF_FDG", "baseline_cn_anticlust")
                opts.fileprefixA {mustBeTextScalar} = "baseline_cn_repA"
                opts.fileprefixB {mustBeTextScalar} = "baseline_cn_repB"
                opts.ARI_thresh {mustBeNumeric} = NaN
                opts.Niter {mustBeNumeric} = 50
                opts.mask mlfourd.ImagingContext2 = ...
                    mlfourd.ImagingContext2( ...
                    fullfile(getenv("SINGULARITY_HOME"), "ADNI", "VolBin", "mask.nii.gz"));
            end

            if isnan(opts.ARI_thresh)
                switch opts.span
                    case 2
                        opts.ARI_thresh = 0.822;
                    case 4
                        opts.ARI_thresh = 0.743;
                    case 6
                        opts.ARI_thresh = 0.78;
                    case 8
                        opts.ARI_thresh = 0.876;
                    case 10
                        opts.ARI_thresh = 0.893;
                    case 12
                        opts.ARI_thresh = 0.903;
                    case 14
                        opts.ARI_thresh = 0.854;
                    case 16
                        opts.ARI_thresh = 0.823;
                    case 18
                        opts.ARI_thresh = 0.836;
                    case 20
                        opts.ARI_thresh = 0.846;
                    case 22
                        opts.ARI_thresh = 0.849;
                    case 24
                        opts.ARI_thresh = 0.889;
                    otherwise
                        error("mladni:ValueError", stackstr())
                end
            end

            % mask_vec = reshape(double(opts.mask), [numel(opts.mask), 1]);
            high_mode_A = [];
            high_mode_B = [];
            low_mode_A = [];
            low_mode_B = [];
            low_mode_indices = [];
            low_mode_ARIs = [];

            for iter = 1:opts.Niter
                try
                    pthA = fullfile(opts.anticlust_home, opts.fileprefixA+iter);
                    pthB = fullfile(opts.anticlust_home, opts.fileprefixB+iter);

                    [A,img_A] = results2num(path2result(pthA));
                    [B,img_B] = results2num(path2result(pthB));

                    ARI = mladni.AdjRandIndex.evaluate_pair(A, B, do_hungarian=false);
                    if ARI > opts.ARI_thresh
                        if isempty(high_mode_A)
                            high_mode_A = img_A;
                        else
                            high_mode_A = cat(5, high_mode_A, img_A);
                        end
                        if isempty(high_mode_B)
                            high_mode_B = img_B;
                        else
                            high_mode_B = cat(5, high_mode_B, img_B);
                        end
                    else
                        if isempty(low_mode_A)
                            low_mode_A = img_A;
                        else
                            low_mode_A = cat(5, low_mode_A, img_A);
                        end
                        if isempty(low_mode_B)
                            low_mode_B = img_B;
                        else
                            low_mode_B = cat(5, low_mode_B, img_B);
                        end

                        low_mode_indices = [low_mode_indices, iter]; %#ok<AGROW>
                        low_mode_ARIs = [low_mode_ARIs, ARI]; %#ok<AGROW>
                    end
                catch ME
                    handwarning(ME)
                end
            end

            %% assign output, S

            S.low_mode_indices = low_mode_indices;
            S.low_mode_ARIs = low_mode_ARIs;

            all = cat(5, high_mode_A, high_mode_B, low_mode_A, low_mode_B);
            S.all_median = stat_of_iter(all, "all_median", @median);
            S.all_iqr = stat_of_iter(all, "all_iqr", @iqr);

            S.high_mode_A_median = stat_of_iter(high_mode_A, "high_mode_A_median", @median);
            S.high_mode_A_iqr = stat_of_iter(high_mode_A, "high_mode_A_iqr", @iqr);
            S.high_mode_B_median = stat_of_iter(high_mode_B, "high_mode_B_median", @median);
            S.high_mode_B_iqr = stat_of_iter(high_mode_B, "high_mode_B_iqr", @iqr);
            S.high_mode_median = stat_of_iter(abs(high_mode_A - high_mode_B), "high_mode_median", @median);
            S.high_mode_iqr = stat_of_iter(abs(high_mode_A - high_mode_B), "high_mode_iqr", @iqr);
            % S.high_wb_median = stat_of_iter(abs(sum(high_mode_A,4) - sum(high_mode_B,4)), "high_wb_median", @median);
            % S.high_wb_iqr = stat_of_iter(abs(sum(high_mode_A,4) - sum(high_mode_B,4)), "high_wb_iqr", @iqr);

            S.low_mode_A_median = stat_of_iter(low_mode_A, "low_mode_A_median", @median);
            S.low_mode_A_iqr = stat_of_iter(low_mode_A, "low_mode_A_iqr", @iqr);
            S.low_mode_B_median = stat_of_iter(low_mode_B, "low_mode_B_median", @median);
            S.low_mode_B_iqr = stat_of_iter(low_mode_B, "low_mode_B_iqr", @iqr);
            S.low_mode_median = stat_of_iter(abs(low_mode_A - low_mode_B), "low_mode_median", @median);
            S.low_mode_iqr = stat_of_iter(abs(low_mode_A - low_mode_B), "low_mode_iqr", @iqr);
            % S.low_wb_median = stat_of_iter(abs(sum(low_mode_A,4) - sum(low_mode_B,4)), "low_wb_median", @median);
            % S.low_wb_iqr = stat_of_iter(abs(sum(low_mode_A,4) - sum(low_mode_B,4)), "low_wb_iqr", @iqr);

            S.delta_delta_of_modes = abs(S.low_mode_median - S.high_mode_median);
            S.delta_delta_of_modes.fileprefix = "delta_delta_of_modes";

            view(S.delta_delta_of_modes, S.all_median, S.all_iqr)

            % high_mode = cat(5, high_mode_A, high_mode_B);
            % low_mode = cat(5, low_mode_A, low_mode_B);
            % S.high_mode_median = stat_of_iter(high_mode, "high_mode_median", @median);
            % S.high_mode_iqr = stat_of_iter(high_mode, "high_mode_iqr", @iqr);
            % S.low_mode_median = stat_of_iter(low_mode, "low_mode_median", @median);
            % S.low_mode_iqr = stat_of_iter(low_mode, "low_mode_iqr", @iqr);
            % 
            % numer = S.low_mode_median - S.high_mode_median;
            % denom = (S.low_mode_median + S.high_mode_median) / 2; 
            % ic = numer ./ denom;
            % ic = ic.scrubNanInf;
            % ic.fileprefix = "residual_mode";
            % view(ic)
            % S.residual_mode = ic;            

            %% inner functions

            function result = path2result(pth)
                mat = fullfile(pth, "NumBases"+opts.span, "OPNMF", "ResultsExtractBases.mat");
                result = load(mat);
            end
            function [B,img] = results2num(results)
                %% sum all patterns from an NMF factoring

                B = results.B;  % Nvoxels x Npatterns
                img = reshape(B, [size(opts.mask), opts.span]);
            end
            function ic = stat_of_iter(img, fp, stat)
                ic = copy(opts.mask);
                ic.selectImagingTool(img=stat(img, 5), fileprefix=fp);
            end

        end
        
        function ri = rand_index(p1, p2, varargin)
            %%  RAND_INDEX Computes the rand index between two partitions.
            %   RAND_INDEX(p1, p2) computes the rand index between partitions p1 and
            %   p2. Both p1 and p2 must be specified as N-by-1 or 1-by-N vectors in
            %   which each elements is an integer indicating which cluster the point
            %   belongs to.
            %
            %   RAND_INDEX(p1, p2, 'adjusted') computes the adjusted rand index
            %   between partitions p1 and p2. The adjustment accounts for chance
            %   correlation.
            %
            %   https://github.com/cmccomb/rand_index/tree/master

            %% Parse the input and throw errors
            % Check inputs
            adj = 0;
            if nargin == 0
                error('Arguments must be supplied.');
            end
            if nargin == 1
                error('Two partitions must be supplied.');
            end
            if nargin > 3
                error('Too many input arguments');
            end
            if nargin == 3
                if strcmp(varargin{1}, 'adjusted')
                    adj = 1;
                else
                    error('%s is an unrecognized argument.', varargin{1});
                end
            end
            if length(p1)~=length(p2)
                error('Both partitions must contain the same number of points.');
            end

            % Check if arguments need to be flattened
            if length(p1)~=numel(p1)
                p1 = p1(:);
                warning('The first partition was flattened to a 1D vector.')
            end
            if length(p2)~=numel(p2)
                p2 = p2(:);
                warning('The second partition was flattened to a 1D vector.')
            end

            % Check for integers
            if isreal(p1) && all(rem(p1, 1)==0)
                % all is well
            else
                warning('The first partition contains non-integers, which may make the results meaningless. Attempting to continue calculations.');
            end

            if isreal(p2) && all(rem(p2, 1)==0)
                % all is well
            else
                warning('The second partition contains non-integers, which may make the results meaningless. Attempting to continue calculations.');
            end

        	%% Preliminary computations and cleansing of the partitions
            N = length(p1);
            [~, ~, p1] = unique(p1);
            N1 = max(p1);
            [~, ~, p2] = unique(p2);
            N2 = max(p2);

            n = zeros(N1,N2);
            %% Create the matching matrix
            for i = 1:length(p1)
                n(p1(i), p2(i)) = n(p1(i), p2(i)) + 1;
            end

            %% If required, calculate the basic rand index
            if adj==0
                ss = sum(sum(n.^2));
                ss1 = sum(sum(n,1).^2);
                ss2 =sum(sum(n,2).^2);
                ri = (nchoosek2(N,2) + ss - 0.5*ss1 - 0.5*ss2)/nchoosek2(N,2);
            end


            %% Otherwise, calculate the adjusted rand index
            if adj==1
                ssm = 0;
                sm1 = 0;
                sm2 = 0;
                for i=1:1:N1
                    for j=1:1:N2
                        ssm = ssm + nchoosek2(n(i,j),2);
                    end
                end
                temp = sum(n,2);
                for i=1:1:N1
                    sm1 = sm1 + nchoosek2(temp(i),2);
                end
                temp = sum(n,1);
                for i=1:1:N2
                    sm2 = sm2 + nchoosek2(temp(i),2);
                end
                NN = ssm - sm1*sm2/nchoosek2(N,2);
                DD = (sm1 + sm2)/2 - sm1*sm2/nchoosek2(N,2);

                % Special case to handle perfect partitions
                if (NN == 0 && DD==0)
                    ri = 1;
                else
                    ri = NN/DD;
                end
            end

            %% Special definition of n choose k
            function c = nchoosek2(a,b)
                if a>1
                    c = nchoosek(a,b);
                else
                    c = 0;
                end
            end
        end

        function test_rand_index()
            %%  https://github.com/cmccomb/rand_index/tree/master

            %% Two arbitrary partitions
            p1 = [1 2 3 3 2 1 1 3 3 1 2 2];
            p2 = [3 2 3 2 2 1 1 2 3 1 3 1];

            % Compute the unadjusted rand index
            ri = mladni.AdjRandIndex.rand_index(p1, p2);
            assert(ri == 7/11);

            % Compute the adjusted rand index
            ari = mladni.AdjRandIndex.rand_index(p1, p2, 'adjusted');
            assert(ari == 1/12);

            %% Two perfect partitions
            p1 = [1 1 1];
            p2 = [1 1 1];

            % Compute the unadjusted rand index
            ri = mladni.AdjRandIndex.rand_index(p1, p2);
            assert(ri == 1);

            % Compute the adjusted rand index
            ari = mladni.AdjRandIndex.rand_index(p1, p2, 'adjusted');
            assert(ari == 1);            
        end

        %% Hungarian algorithm/method.
        %  See also https://en.wikipedia.org/wiki/Hungarian_algorithm

        function [Matching,Cost] = Hungarian(Perf)
            %
            % [MATCHING,COST] = Hungarian_New(WEIGHTS)
            %
            % A function for finding a minimum edge weight matching given a MxN Edge
            % weight matrix WEIGHTS using the Hungarian Algorithm.
            %
            % An edge weight of Inf indicates that the pair of vertices given by its
            % position have no adjacent edge.
            %
            % MATCHING return a MxN matrix with ones in the place of the matchings and
            % zeros elsewhere.
            %
            % COST returns the cost of the minimum matching

            % Written by: Alex Melin 30 June 2006

            % Initialize Variables
            Matching = zeros(size(Perf));

            % Condense the Performance Matrix by removing any unconnected vertices to
            % increase the speed of the algorithm

            % Find the number in each column that are connected
            num_y = sum(~isinf(Perf),1);
            % Find the number in each row that are connected
            num_x = sum(~isinf(Perf),2);

            % Find the columns(vertices) and rows(vertices) that are isolated
            x_con = find(num_x~=0);
            y_con = find(num_y~=0);

            % Assemble Condensed Performance Matrix
            P_size = max(length(x_con),length(y_con));
            P_cond = zeros(P_size);
            P_cond(1:length(x_con),1:length(y_con)) = Perf(x_con,y_con);
            if isempty(P_cond)
                Cost = 0;
                return
            end

            % Ensure that a perfect matching exists
            % Calculate a form of the Edge Matrix
            Edge = P_cond;
            Edge(P_cond~=Inf) = 0;
            % Find the deficiency(CNUM) in the Edge Matrix
            cnum = mladni.AdjRandIndex.min_line_cover(Edge);

            % Project additional vertices and edges so that a perfect matching
            % exists
            Pmax = max(max(P_cond(P_cond~=Inf)));
            P_size = length(P_cond)+cnum;
            P_cond = ones(P_size)*Pmax;
            P_cond(1:length(x_con),1:length(y_con)) = Perf(x_con,y_con);

            %*************************************************
            % MAIN PROGRAM: CONTROLS WHICH STEP IS EXECUTED
            %*************************************************
            exit_flag = 1;
            stepnum = 1;
            while exit_flag
                switch stepnum
                    case 1
                        [P_cond,stepnum] = mladni.AdjRandIndex.step1(P_cond);
                    case 2
                        [r_cov,c_cov,M,stepnum] = mladni.AdjRandIndex.step2(P_cond);
                    case 3
                        [c_cov,stepnum] = mladni.AdjRandIndex.step3(M,P_size);
                    case 4
                        [M,r_cov,c_cov,Z_r,Z_c,stepnum] = mladni.AdjRandIndex.step4(P_cond,r_cov,c_cov,M);
                    case 5
                        [M,r_cov,c_cov,stepnum] = mladni.AdjRandIndex.step5(M,Z_r,Z_c,r_cov,c_cov);
                    case 6
                        [P_cond,stepnum] = mladni.AdjRandIndex.step6(P_cond,r_cov,c_cov);
                    case 7
                        exit_flag = 0;
                end
            end

            % Remove all the virtual satellites and targets and uncondense the
            % Matching to the size of the original performance matrix.
            Matching(x_con,y_con) = M(1:length(x_con),1:length(y_con));
            Cost = sum(sum(Perf(Matching==1)));
        end
        function [P_cond,stepnum] = step1(P_cond)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   STEP 1: Find the smallest number of zeros in each row
            %           and subtract that minimum from its row
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            P_size = length(P_cond);

            % Loop throught each row
            for ii = 1:P_size
                rmin = min(P_cond(ii,:));
                P_cond(ii,:) = P_cond(ii,:)-rmin;
            end

            stepnum = 2;
        end
        function [r_cov,c_cov,M,stepnum] = step2(P_cond)
            %**************************************************************************
            %   STEP 2: Find a zero in P_cond. If there are no starred zeros in its
            %           column or row start the zero. Repeat for each zero
            %**************************************************************************

            % Define variables
            P_size = length(P_cond);
            r_cov = zeros(P_size,1);  % A vector that shows if a row is covered
            c_cov = zeros(P_size,1);  % A vector that shows if a column is covered
            M = zeros(P_size);        % A mask that shows if a position is starred or primed

            for ii = 1:P_size
                for jj = 1:P_size
                    if P_cond(ii,jj) == 0 && r_cov(ii) == 0 && c_cov(jj) == 0
                        M(ii,jj) = 1;
                        r_cov(ii) = 1;
                        c_cov(jj) = 1;
                    end
                end
            end

            % Re-initialize the cover vectors
            r_cov = zeros(P_size,1);  % A vector that shows if a row is covered
            c_cov = zeros(P_size,1);  % A vector that shows if a column is covered
            stepnum = 3;
        end
        function [c_cov,stepnum] = step3(M,P_size)
            %**************************************************************************
            %   STEP 3: Cover each column with a starred zero. If all the columns are
            %           covered then the matching is maximum
            %**************************************************************************

            c_cov = sum(M,1);
            if sum(c_cov) == P_size
                stepnum = 7;
            else
                stepnum = 4;
            end
        end
        function [M,r_cov,c_cov,Z_r,Z_c,stepnum] = step4(P_cond,r_cov,c_cov,M)
            %**************************************************************************
            %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
            %           zero in the row containing this primed zero, Go to Step 5.
            %           Otherwise, cover this row and uncover the column containing
            %           the starred zero. Continue in this manner until there are no
            %           uncovered zeros left. Save the smallest uncovered value and
            %           Go to Step 6.
            %**************************************************************************

            P_size = length(P_cond);

            zflag = 1;
            while zflag
                % Find the first uncovered zero
                row = 0; col = 0; exit_flag = 1;
                ii = 1; jj = 1;
                while exit_flag
                    if P_cond(ii,jj) == 0 && r_cov(ii) == 0 && c_cov(jj) == 0
                        row = ii;
                        col = jj;
                        exit_flag = 0;
                    end
                    jj = jj + 1;
                    if jj > P_size; jj = 1; ii = ii+1; end
                    if ii > P_size; exit_flag = 0; end
                end

                % If there are no uncovered zeros go to step 6
                if row == 0
                    stepnum = 6;
                    zflag = 0;
                    Z_r = 0;
                    Z_c = 0;
                else
                    % Prime the uncovered zero
                    M(row,col) = 2;
                    % If there is a starred zero in that row
                    % Cover the row and uncover the column containing the zero
                    if sum(find(M(row,:)==1)) ~= 0
                        r_cov(row) = 1;
                        zcol = find(M(row,:)==1);
                        c_cov(zcol) = 0;
                    else
                        stepnum = 5;
                        zflag = 0;
                        Z_r = row;
                        Z_c = col;
                    end
                end
            end
        end
        function [M,r_cov,c_cov,stepnum] = step5(M,Z_r,Z_c,r_cov,c_cov)
            %**************************************************************************
            % STEP 5: Construct a series of alternating primed and starred zeros as
            %         follows.  Let Z0 represent the uncovered primed zero found in Step 4.
            %         Let Z1 denote the starred zero in the column of Z0 (if any).
            %         Let Z2 denote the primed zero in the row of Z1 (there will always
            %         be one).  Continue until the series terminates at a primed zero
            %         that has no starred zero in its column.  Unstar each starred
            %         zero of the series, star each primed zero of the series, erase
            %         all primes and uncover every line in the matrix.  Return to Step 3.
            %**************************************************************************

            zflag = 1;
            ii = 1;
            while zflag
                % Find the index number of the starred zero in the column
                rindex = find(M(:,Z_c(ii))==1);
                if rindex > 0
                    % Save the starred zero
                    ii = ii+1;
                    % Save the row of the starred zero
                    Z_r(ii,1) = rindex;
                    % The column of the starred zero is the same as the column of the
                    % primed zero
                    Z_c(ii,1) = Z_c(ii-1);
                else
                    zflag = 0;
                end

                % Continue if there is a starred zero in the column of the primed zero
                if zflag == 1
                    % Find the column of the primed zero in the last starred zeros row
                    cindex = find(M(Z_r(ii),:)==2);
                    ii = ii+1;
                    Z_r(ii,1) = Z_r(ii-1);
                    Z_c(ii,1) = cindex;
                end
            end

            % UNSTAR all the starred zeros in the path and STAR all primed zeros
            for ii = 1:length(Z_r)
                if M(Z_r(ii),Z_c(ii)) == 1
                    M(Z_r(ii),Z_c(ii)) = 0;
                else
                    M(Z_r(ii),Z_c(ii)) = 1;
                end
            end

            % Clear the covers
            r_cov = r_cov.*0;
            c_cov = c_cov.*0;

            % Remove all the primes
            M(M==2) = 0;

            stepnum = 3;

        end
        function [P_cond,stepnum] = step6(P_cond,r_cov,c_cov)
            % *************************************************************************
            % STEP 6: Add the minimum uncovered value to every element of each covered
            %         row, and subtract it from every element of each uncovered column.
            %         Return to Step 4 without altering any stars, primes, or covered lines.
            %**************************************************************************

            a = find(r_cov == 0);
            b = find(c_cov == 0);
            minval = min(min(P_cond(a,b)));

            P_cond(find(r_cov == 1),:) = P_cond(find(r_cov == 1),:) + minval;
            P_cond(:,find(c_cov == 0)) = P_cond(:,find(c_cov == 0)) - minval;

            stepnum = 4;
        end
        function cnum = min_line_cover(Edge)
            % Step 2
            [r_cov,c_cov,M,stepnum] = mladni.AdjRandIndex.step2(Edge);
            % Step 3
            [c_cov,stepnum] = mladni.AdjRandIndex.step3(M,length(Edge));
            % Step 4
            [M,r_cov,c_cov,Z_r,Z_c,stepnum] = mladni.AdjRandIndex.step4(Edge,r_cov,c_cov,M);
            % Calculate the deficiency
            cnum = length(Edge)-sum(r_cov)-sum(c_cov);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end