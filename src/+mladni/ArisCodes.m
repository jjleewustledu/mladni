classdef ArisCodes < handle
    %% https://github.com/sotiraslab/aris_nmf_analyses
    %  
    %  Created 20-Apr-2024 13:37:17 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2568132 (R2024a) Update 1 for MACA64.  Copyright 2024 John J. Lee.
    
    methods (Static)
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
            % opts.K: count of model spans, e.g., K == 20 examines model spans 2:2:40
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
                opts.K {mustBeInteger} = 20
                opts.do_plot logical = true
            end
            pathDir1 = convertStringsToChars(pathDir1);
            pathDir2 = convertStringsToChars(pathDir2);
            outputDir = convertStringsToChars(outputDir);

            K=opts.K;
            
            sortedBasisNum=2:2:2*K;
            
            ARI=zeros(K,1);
            
            for exp_=1:K
                disp([ num2str(exp_) '/' num2str(K)])

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
                    ARI(exp_) = mladni.NMF.clustering_adjustedRand_fast(clustering1,clustering2);                
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
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
