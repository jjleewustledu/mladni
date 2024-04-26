classdef ArisCodes < handle
    %% https://github.com/sotiraslab/aris_nmf_analyses
    %  
    %  Created 20-Apr-2024 13:37:17 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2568132 (R2024a) Update 1 for MACA64.  Copyright 2024 John J. Lee.
    
    methods (Static)
        function [RecError,X] = calcBLSARecError(inFiles, resultDir)

            resultDir = convertStringsToChars(resultDir);

            % populating path
            % addpath(genpath('../../'));

            mfMeth = {'OPNMF','PCA','Jade'} ;

            % load data -> we need to estimate the age
            param.isList = 1 ;
            param.downSample = 1 ;
            data = mladni.ArisCodes.loadData(inFiles,param,[]) ;
            meanX = mean(data.X,2);
            X = data.X(meanX>0,:); clear data

            % hard coded results path
            listing = dir(resultDir);
            listing=listing(3:end) ;
            hh =cellfun(@(x) (strfind(x,'NumBases')==1),{listing(:).name},'UniformOutput',false) ;
            listing=listing(cellfun(@(x) (~isempty(x)&&(x==1)),hh));
            numDifBases=numel(listing) ;

            % sort them in ascending order
            basisNum = zeros(1,numDifBases) ;
            for i=1:numDifBases
                basisNum(i) = str2double(listing(i).name(9:end));    
            end
            [~,idx]=sort(basisNum) ;
            sortedBasisNum=basisNum(idx) ;

            % allocate memory 
            RecError=zeros(numDifBases,numel(mfMeth));

            for j=1:numel(mfMeth)
                for i=1:numDifBases
                    dataPath=[resultDir '/' listing(idx(i)).name '/' mfMeth{j} '/ResultsExtractBases.mat'];
                    % loading B (bases) and C (loading coefficients) - we do not
                    % need B
                    load(dataPath); %#ok<LOAD>
                    if(strcmp(mfMeth{j},'PCA'))
                        Est = B(meanX>0,:)*C + meanX(meanX>0)*ones(1,size(C,2)); clear B C
                    else
                        Est = B(meanX>0,:)*C ; clear B C
                    end
                    RecError(i,j) = norm(X-Est,'fro') ;
                    clear Est
                end
            end            
        end

        function RecError = calcRecError(X, resultDir, resultsDir, selectVoxels, opts)
            
            % Function that calculates the reconstruction error given data matrix X and the non-negative matrix 
            % factorization results saved in the resultDir directory. The function additionally saves figures that plot
            % the reconstruction error, the gradient of the reconstruction error, and the percentage of improvement as a 
            % function of the number of components
            
            % load results and calculate reconstruction error
            % We assume that the results directory is organized as follows: there is a folder for every solution (i.e.,
            % NumBases${K}, wehere K is the number of components for the solution), and inside the folder there is .mat 
            % file (ResultsExtractBases.mat) that contains the matrices W and H that were estimated by the non negative 
            % matrix factorization

            % selectVoxels is a logical vector of size N_voxels x 1, for N_voxels in the imaging field of view

            arguments
                X {mustBeNumeric}
                resultDir {mustBeFolder}
                resultsDir {mustBeTextScalar}
                selectVoxels logical
                opts.do_plot logical = false
            end

            resultDir = convertStringsToChars(resultDir); %% JJL
            resultsDir = convertStringsToChars(resultsDir); %% JJL    
            assert(islogical(selectVoxels))
            assert(size(X,1) == sum(selectVoxels, "all"), stackstr()) %% JJL's BUG CHECK
            
            listing = dir(resultDir);
            listing=listing(3:end);
            hh =cellfun(@(x) ( (strfind(x,'NumBases')==1)  ),{listing(:).name},'UniformOutput',false) ;
            listing=listing(cellfun(@(x) (~isempty(x)&&(x==1)),hh));
            numDifBases=numel(listing);
            
            % sort them in ascending order
            basisNum = zeros(1,numDifBases) ;
            for i=1:numDifBases
                basisNum(i) = str2double(listing(i).name(9:end));
            end
            [~,idx]=sort(basisNum) ;
            sortedBasisNum=basisNum(idx) ;
            
            RecError=zeros(numDifBases,1);
            for b=1:numDifBases
                disp(b/numDifBases)
                load( ...
                    fullfile(resultDir, ['NumBases', num2str(sortedBasisNum(b))], 'OPNMF', ...
                    'ResultsExtractBases.mat')) %#ok<LOAD>
                Est = B*C ;
                assert(size(Est,1) == size(selectVoxels,1), stackstr()) %% JJL's BUG CHECK
                Est = Est(selectVoxels,:); %% JJL
                RecError(b) = norm(X-Est,'fro') ;
                clear B C
            end
            
            if opts.do_plot
                plt = mladni.Plot();

                % make figures
                % 1) reconstruction error
                plt.plotxy(sortedBasisNum, RecError, ...
                    xlab="Number of patterns in model space", ...
                    ylab="Reconstruction error", ...
                    fileprefix=fullfile(resultsDir, "RecError"));

                % 2) gradient of reconstruction error
                plt.plotxy(sortedBasisNum(2:end), diff(RecError), ...
                    xlab="Number of patterns in model space", ...
                    ylab="Gradient of reconstruction error", ...
                    fileprefix=fullfile(resultsDir, "gradientRecError"));

                % 3) Fractional improvement over range of components used
                improvement = abs(RecError-RecError(1))./abs(RecError(1)-RecError(end));
                plt.plotxy(sortedBasisNum, improvement, ...
                    xlab="Number of patterns in model space", ...
                    ylab="Improvement with more patterns", ...
                    fileprefix=fullfile(resultsDir, "fractionalImprovementRecError"));

                %close all
            end
        end        
        
        function calculateSelectedComponentWeightedAverageNIFTI(dataLists,basesDir,numBases,outFile)
            %% calculates projections of NIfTI to each of the bases in basesDir.
            %  Selects single folder NumBases\d, then accesses all /OPNMF/niiImg/Basis_*.nii to project files in dataList to Basis_*.nii
            %
            % dataList: .csv file containing the images for which the coefficients need
            %           to be calculated. Full path for every image file is given in
            %           every line
            % resultsDir: directory where the NMF results are (i.e., the level where
            %             the NumBases folder is placed)
            % numBases: determines the solution for which one wants to calculate
            %           subject coefficients
            % outFile: CSV-file saving output
            
            import mlniftitools.*;
            dataLists = convertStringsToChars(dataLists);
            basesDir = convertStringsToChars(basesDir);
            outFile = convertStringsToChars(outFile);
        
            % works with only single folder NumBases\d found from reproducibility analysis

            fprintf("%s: numBases->%g\n", stackstr(), numBases)
            dataPath=[basesDir '/NumBases' num2str(numBases) '/OPNMF/niiImg/'];
            
            % loading estimated non-negative components
            listing = glob(fullfile(dataPath, 'Basis_*.nii'));
            listing = listing(~contains(listing, 'argmax') & ~contains(listing, 'all'));
            
            if(length(listing)~=numBases)
                error(['I cannot find ' num2str(numBases) ' basis images in the folder ' dataPath ])
            end
            
            for i=1:numBases
                nii = load_untouch_nii(listing{i});
                B(:,i) = double(nii.img(:)');
            end
            
            % normalize to sum to 1
            Blen = sum(B,1);
            if any(Blen==0)
                Blen(Blen==0) = 1;
            end
            nB = bsxfun(@times,B,1./Blen);
            
            % since the size and number of files is such that we can not manage
            % in batch mode, we are going to calculate weighted average values
            % subject by subject
            dataLists = ensureCell(dataLists);
            for dl = 1:length(dataLists)
    
                % read list
                fid=fopen(dataLists{dl},'r');
                if (fid == -1)
                    error(['extractBases:calculateComponentWeightedAverage ','Can not open ' list ' file.']);
                end
                datafullpath = textscan(fid,'%s\n');
                fclose(fid);
                
                datafullpath = datafullpath{1,1} ;
                datafullpath = cellstr(datafullpath) ;
                count = numel(datafullpath);
                
                fid = fopen(outFile,'w');
                frmtWrite='%s,';
                frmtWrite=[frmtWrite repmat('%f,',1,numBases-1)]; %#ok<*AGROW>
                frmtWrite=[frmtWrite '%f\n'];
                
                wA = zeros(count,numBases) ;

                for i=1:count                    
                    %print statement suggested by tom
                    if ~mod(i,10)
                        disp(sprintf('%s: (%d/%d)', stackstr(), i, count)); %#ok<DSPS>
                    end
                    nii = load_untouch_nii(datafullpath{i});
                    wA(i,:) = double(nii.img(:)')*nB;
                    fprintf(fid,frmtWrite,datafullpath{i},wA(i,:)');
                end
                fclose(fid);
            end
        end

        function r = clustering_adjustedRand_fast(u,v)
            % clustering quality measures assumptions : 
            % 0 corresponds to background
            % we do not care about the background
            % cluster labels are assumed to be enumerated from 1 to max number of
            % clusters
            %
            % this function should not be used when comparing binary segmentations
            
            m=max(max(u),max(v));
            
            if(m == 1)
                error('ClusteringAdjustedRandFast:argChk','This method should not be used for comparing binary segmentations');
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

        function [meanInner,medianInner,ARI,overlap,sortedBasisNum] = evaluateReproducibility( ...
                pathDir1, pathDir2, outputDir, opts)
            %% https://github.com/sotiraslab/aris_nmf_analyses/blob/main/evaluateReproducibility.m
            %
            % This function evaluates the reproducibility between the results obtained
            % for two parts of a data split.
            % The scripts assumes that experiments using the same range of components
            % have been performed for both splits. 
            %
            % inputs
            % pathDir1: path to the results of the first split
            % pathDir2: path to the results of the second split
            % outputDir: where to save results and figures
            % opts.N_K: count of models with incremented cardinality, e.g., N_K == 20 examines model spans 2:2:40
            % opts.do_plot: do plot & save figures
            %
            % outputs
            % meanInner : mean value of the inner product between matched components
            % medianInner : median value of the inner product between matched components
            % ARI : adjusted Rand Index evaluated by deriving hard clusters from the estimated components
            % overlap : double, size ~ {1,numDifBases}[length(wlen1) in 2:2:40, 1]
            % sortedBasisNum : size ~ [numDifBases,1]
            
            arguments
                pathDir1 char {mustBeTextScalar,mustBeFolder}
                pathDir2 char {mustBeTextScalar,mustBeFolder}
                outputDir char {mustBeTextScalar}
                opts.N_K {mustBeInteger} = 20
                opts.do_plot logical = true
                opts.rep double = 1
                opts.N_rep double = 1
            end
            pathDir1 = convertStringsToChars(pathDir1);
            pathDir2 = convertStringsToChars(pathDir2);
            outputDir = convertStringsToChars(outputDir);
            
            sortedBasisNum=2:2:2*opts.N_K;
            
            ARI=zeros(opts.N_K,1);
            
            for exp_=1:opts.N_K
                try               
                    resSplit1 = load([pathDir1 '/NumBases' num2str(sortedBasisNum(exp_)) '/OPNMF/ResultsExtractBases.mat']) ;
                    resSplit2 = load([pathDir2 '/NumBases' num2str(sortedBasisNum(exp_)) '/OPNMF/ResultsExtractBases.mat']) ;

                    % normalize to unit norm
                    wlen1 = sqrt(sum((resSplit1.B).^2)) ;  % B ~ Nvoxels x Ncomponents
                    wlen2 = sqrt(sum((resSplit2.B).^2)) ;  % wlen{1,2} ~ 1 x Ncomponents

                    if any(wlen1==0)
                        wlen1(wlen1==0) = 1;
                    end
                    W1 = bsxfun(@times,resSplit1.B,1./wlen1) ;

                    if any(wlen2==0)
                        wlen2(wlen2==0) = 1;
                    end
                    W2 = bsxfun(@times,resSplit2.B,1./wlen2) ;

                    % calculate inner products
                    inner_product = W1'*W2 ;  % Ncomponents x Ncomponents

                    % take a distance
                    dist = 2*(1 - inner_product) ;

                    % find correspondences
                    [Matching,~] = Hungarian(dist);
                    [~,idx_hug1]=max(Matching,[],2);

                    % overlap - hungarian
                    overlap{exp_} = zeros(length(wlen1),1) ;
                    for b=1:length(wlen1)
                        overlap{exp_}(b) = inner_product(b,idx_hug1(b));
                    end

                    % overlap with best
                    overlap_best{exp_} = max(inner_product,[],2) ;

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
                    ARI(exp_) = mladni.ArisCodes.clustering_adjustedRand_fast(clustering1,clustering2);     

                    fexp_ = exp_ / opts.N_K;
                    frep = ((opts.rep - 1) + fexp_) / opts.N_rep;
                    progressbar(frep, fexp_)
                catch ME
                    handwarning(ME);
                end
            end
            
            meanInner=cell2mat(cellfun(@(x) mean(x),overlap,'UniformOutput', false));
            medianInner=cell2mat(cellfun(@(x) median(x),overlap,'UniformOutput', false));
            
            % save results
            save([outputDir '/reproducibilityResults.mat'],'meanInner','medianInner','ARI','sortedBasisNum');
            if opts.do_plot
                %stdInner=cellfun(@(x) std(x),overlap,'UniformOutput', false);
                figure;plot(sortedBasisNum,meanInner,'b','LineWidth',2)
                xlabel('Number of components','fontsize',12)
                ylabel({'Split-sample reproducibility';'(mean inner product)'},'fontsize',12)
                set(gca,'fontsize',12)
                saveas(gcf,[outputDir '/MeanInnerProductReproducibility.fig'])
                saveas(gcf,[outputDir '/MeanInnerProductReproducibility.png'])
                
                figure;plot(sortedBasisNum,medianInner,'b','LineWidth',2)
                xlabel('Number of components','fontsize',12)
                ylabel({'Split-sample reproducibility';'(median inner product)'},'fontsize',12)
                set(gca,'fontsize',12)
                saveas(gcf,[outputDir '/MedianInnerProductReproducibility.fig'])
                saveas(gcf,[outputDir '/MedianInnerProductReproducibility.png'])
                
                figure;plot(sortedBasisNum,ARI,'b','LineWidth',2)
                xlabel('Number of components','fontsize',12)
                ylabel({'Split-sample reproducibility' ;'(Adjusted Rand Index)'},'fontsize',12)
                set(gca,'fontsize',12)
                saveas(gcf,[outputDir '/AdjustedRandIndexReproducibility.fig'])
                saveas(gcf,[outputDir '/AdjustedRandIndexReproducibility.png'])
            end
        end

        function data = loadData(dataInput,param,subsetIdx)
            
            import mlniftitools.*;

            dataInput = convertStringsToChars(dataInput);
            
            if(param.isList==1)
                disp(['List ' dataInput ]);
                data = mladni.ArisCodes.loadDataFromList(dataInput,param,subsetIdx) ; % returned data is a struct !
            else
                % if data input is not given as list, then what we do depends on
                % whether the code is used as a function or a standalone executable
                if(isdeployed)
                    % dataInput contains the full path to the .mat that contains
                    % data
                    tmp = load(dataInput,'data');
                    tmp = struct2cell(tmp);
                    data = tmp{1}; clear tmp
                    data.dsflag = 0;
                else
                    % dataInput is the data_matrix
                    data = dataInput ; clear dataInput
                end
                % we have loaded data but we don't know if it is a struct or just an
                % array
                if (isstruct(data))
                    % check if field X is there
                    if(~isfield(data,'X'))
                        error('extractBases:argChk','The field X that should contain the data is missing !');
                    end
                    if(~isempty(subsetIdx))
                        data.X = data.X(:,subsetIdx==1);
                    end
                    if(isfield(data,'dimx') && isfield(data,'dimy') && isfield(data,'dimz') && (param.downSample~=1))
                        % downsample data
                        tmp = data.X ;
                        h_nsamples = size(data.X,2) ;
                        data.X = [] ;
                        new_size_x = ceil(data.dimx/param.downSample) ;
                        new_size_y = ceil(data.dimy/param.downSample) ;
                        new_size_z = ceil(data.dimz/param.downSample) ;
                        data.X = zeros(new_size_x*new_size_y*new_size_z,h_nsamples) ;
                        for i=1:h_nsamples
                            xx1 = linspace(1,data.dimx,data.dimx) ;
                            yy1 = linspace(1,data.dimy,data.dimy) ;
                            zz1 = linspace(1,data.dimz,data.dimz) ;
                            xx2 = linspace(1,data.dimx,new_size_x) ;
                            yy2 = linspace(1,data.dimy,new_size_y) ;
                            zz2 = linspace(1,data.dimz,new_size_z) ;
                            [XX1,YY1,ZZ1] = meshgrid(yy1,xx1,zz1) ;
                            [XX2,YY2,ZZ2] = meshgrid(yy2,xx2,zz2) ;
                            h_img = reshape(tmp(:,i),data.dimx.data.dimy,data.dimz) ;
                            img_interp = interp3(XX1,YY1,ZZ1,h_img,XX2,YY2,ZZ2) ;
                            data.X(:,i) = double(img_interp(:));
                        end
                        clear tmp
                        clear h_nsamples
                        clear h_img
                        data.dsflag = 1;
                    else
                        data.dsflag = 0;
                    end
                else
                    % put data into a struct
                    tmp = data ; clear data ; data.X = tmp ; clear tmp ;
                    if(~isempty(subsetIdx))
                        data.X = data.X(:,subsetIdx==1);
                    end
                    data.dsflag = 0;
                end
            end
            
            if(param.permute)
                s = RandStream.create('mt19937ar','seed',param.permSeed);RandStream.setGlobalStream(s);
                for i=1:size(data.X,2)
                    tmp = data.X(randperm(length(data.X(:,i))),i);
                    data.X(:,i) = tmp ;
                end
            end
        end

        function data = loadDataFromList(list,param,subsetIdx)
            
            % Loading data form list
            % Input:
            %   list      : .txt file that contains function the full path for the data
            %               to be read
            %   param     : structure containing parameters
            %   subsetIdx : either EMPTY or FILENAME of the .txt file that contains
            %               that information regarding which images should be read
            %
            % This function returns a structure array data that has the following fields
            %   X         : matlab array that contains the actual data vectorized
            %               (dimension is D x N -> dimensionality x number of samples)
            %   y         : label information
            %   dimx      : dimension in x-axis
            %   dimy      : dimension in y-axis
            %   dimz      : dimension in z-axis
            %   nii       : nii structure of the input data that may be used in case
            %               one want to save the files in .nii format
            
            import mlniftitools.*;
            
            list = convertStringsToChars(list);
            fid=fopen(list,'r');
            if (fid == -1)
                error( ...
                    "mladni:IOError", ...
                    "extractBases:loadDataFromList cannot open %s.", list);
            end
            datafullpath = textscan(fid,'%s %d\n');
            fclose(fid);
            
            data.y = datafullpath{1,2} ;
            datafullpath = datafullpath{1,1} ;
            datafullpath = cellstr(datafullpath) ;
            
            % keep filename info in the param structure, so that it can be saved
            data.filenames = datafullpath ;
            
            % find indexed entries, if trainIdx is not empty
            if(~isempty(subsetIdx))
                datafullpath = datafullpath(subsetIdx==1) ;
            end
            
            count = numel(datafullpath);
            info = load_untouch_header_only(datafullpath{1});
            
            data.dimx = info.dime.dim(2) ;
            data.dimy = info.dime.dim(3) ;
            data.dimz = info.dime.dim(4) ;
            
            % initialize data structure (vectorize data)
            if(param.downSample==1)
                % if no downsampling is performed
                disp('No downsampling is going to be performed to the data');
                % if data need to be masked, do not initialize yet the data structure
                if(isempty(param.mask))
                    data.X = zeros(data.dimx*data.dimy*data.dimz,count);
                end
                data.dsflag = 0 ;
            else
                % case of downsampling
                disp(['Data are going to be downsampled by a factor of ' num2str(param.downSample)]);
                new_size_x = ceil(data.dimx/param.downSample) ;
                new_size_y = ceil(data.dimy/param.downSample) ;
                new_size_z = ceil(data.dimz/param.downSample) ;
                % if data need to be masked, do not initialize yet the data structure
                if(isempty(param.mask))
                    data.X = zeros(new_size_x*new_size_y*new_size_z,count);
                end
                data.dsflag = 1;
            end
            
            % if data need to smoothed
            if(~isempty(param.smooth) && param.smooth ~= 0)
                disp('Data are going to be smoothed !');
                disp(['User specified smoothing kernel : ' num2str(param.smooth)]);
                [status,scratchDir]=tempdir; %system('echo ${SBIA_TMPDIR}');
                scratchDir=strtrim(scratchDir);
                if(status ~= 0)
                    error(['extractBases:loscratchDiradDataFromList ','Can not find which is the scratch directory.']);
                end
                if(isempty(scratchDir))
                    scratchDir='/tmp';
                end
                [status,pid]=system('echo $$');
                pid=str2double(pid);
                if(status ~= 0)
                    error( ...
                        "mladni:RuntimeError", ...
                        "extractBases:loadDataFromList cannot find process id.");
                end
                % make temp directory to smooth images
                tmpDirName=[scratchDir '/smoothImg_NumBases_' num2str(param.numBase) '_PID_' num2str(pid) ] ;
                
                for tt=1:3
                    % if directory exists I will try to create a unique name 3 times
                    if(exist(tmpDirName,'dir'))
                        tmpDirName = [tmpDirName '_' num2str(floor(rand(1)*10000)) ] ;
                    end
                end
                if(~exist(tmpDirName,'dir'))
                    status=mkdir(tmpDirName);
                    if(status ~= 1)
                        error( ...
                            "mladni:FilesystemError", ...
                            "extractBases:loadDataFromList cannot create temporary directory %s to store smoothed images.", ...
                            tmpDirName);
                    end
                else
                    error("mladni:IOError", ...
                        "extractBases:loadDataFromList:  someone else is writing in temporary directory %s.", ...
                        tmpDirName);
                end
            else
                disp('No additional user specified smoothing is going to performed to the data');
            end
            
            % if data need to be masked
            if(~isempty(param.mask))
                disp('The data are going to be masked');
                disp(['Mask provided by user in ' param.mask]);
                mask_nii = load_untouch_nii(param.mask);
                mask = mask_nii.img ;
                
                % check that mask is given in the correct space
                if(mask_nii.hdr.dime.dim(2) ~= data.dimx || ...
                        mask_nii.hdr.dime.dim(3) ~= data.dimy || ...
                        mask_nii.hdr.dime.dim(4) ~= data.dimz )
                    error("mladni:ValueError", ...
                        "extractBases:loadDataFromList: mask size does not agree with image size.");
                end
                
                % initialize data structure
                if(param.downSample==1)
                    data.X = zeros(sum(mask(:)>0),count);
                else
                    xx1 = linspace(1,data.dimx,data.dimx) ;
                    yy1 = linspace(1,data.dimy,data.dimy) ;
                    zz1 = linspace(1,data.dimz,data.dimz) ;
                    xx2 = linspace(1,data.dimx,new_size_x) ;
                    yy2 = linspace(1,data.dimy,new_size_y) ;
                    zz2 = linspace(1,data.dimz,new_size_z) ;
                    [XX1,YY1,ZZ1] = meshgrid(yy1,xx1,zz1) ;
                    [XX2,YY2,ZZ2] = meshgrid(yy2,xx2,zz2) ;
                    mask = interp3(XX1,YY1,ZZ1,double(mask),XX2,YY2,ZZ2) ;
                    data.X = zeros(sum(mask(:)>0),count) ;
                end
            else
                disp('No user provided mask is specified.');
                mask = [] ;
            end
            data.mask = mask ;
            
            % load data
            % added print statement - Tom
            disp(sprintf("Loading %d images from %s ...", count, list)); %#ok<DSPS>
            for i=1:count                          
                try
                    if(~isempty(param.smooth) && param.smooth ~= 0)
                        smooth_command = ...
                            ['3dmerge -1blur_fwhm ' num2str(param.smooth) ' -prefix ' tmpDirName '/smoothed_image_' num2str(i) '.nii.gz ' datafullpath{i} ] ;
                        [status,~] = system(smooth_command) ;
                        if(status ~= 0)
                            error("mladni:RuntimeError", ...
                                "extractBases:loadDataFromList cannot smooth image %s.", datafullpath{i});
                        end
                        % load smoothed image
                        nii = load_untouch_nii([tmpDirName '/smoothed_image_' num2str(i) '.nii.gz']);

                        % remove smoothed image
                        delete([tmpDirName '/smoothed_image_' num2str(i) '.nii.gz']);
                    else
                        nii = load_untouch_nii(datafullpath{i});
                    end
                    % if downsampling should be performed
                    if(param.downSample==1)
                        if(~isempty(mask))
                            data.X(:,i) = nii.img(mask>0) ;
                        else
                            data.X(:,i) = nii.img(:) ;
                        end
                    else
                        xx1 = linspace(1,size(nii.img,1),size(nii.img,1)) ;
                        yy1 = linspace(1,size(nii.img,2),size(nii.img,2)) ;
                        zz1 = linspace(1,size(nii.img,3),size(nii.img,3)) ;
                        xx2 = linspace(1,size(nii.img,1),new_size_x) ;
                        yy2 = linspace(1,size(nii.img,2),new_size_y) ;
                        zz2 = linspace(1,size(nii.img,3),new_size_z) ;
                        [XX1,YY1,ZZ1] = meshgrid(yy1,xx1,zz1) ;
                        [XX2,YY2,ZZ2] = meshgrid(yy2,xx2,zz2) ;
                        img_interp = interp3(XX1,YY1,ZZ1,double(nii.img),XX2,YY2,ZZ2) ;
                        if(~isempty(mask))
                            h_nii = double(img_interp(mask>0));
                            data.X(:,i) = h_nii(:) ; clear h_nii ;
                        else
                            data.X(:,i) = double(img_interp(:));
                        end
                    end
                catch ME
                    handwarning(ME)
                end
            end
            data.nii = nii ;
            
            % % if data were smoothed, remove directory
            if(~isempty(param.smooth) && param.smooth ~= 0)
                [status,~,~]=rmdir(tmpDirName);
                if(status ~= 1)
                    error("mladni:FilesystemError", ...
                        "extractBases:loadDataFromList cannot remove folder %s.", tmpDirName);
                end
            end
        end
        
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
