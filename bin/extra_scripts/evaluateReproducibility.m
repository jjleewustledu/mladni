function [meanInner,medianInner,ARI,sortedBasisNum] = evaluateReproducibility(pathDir1,pathDir2,outputDir)

% function that evaluates the reproducibility between the results obtained
% for two different splits
% The scripts assumes that experiments using the same range of components
% have been performed for both splits. 

% inputs
% pathDir1: path to the results of the first split
% pathDir2: path to the results of the second split
% outputDir: where to save results and figures

% outputs
% meanInner : mean value of the inner product between matched components
% medianInner : median value of the inner product between matched components
% ARI : adjusted Rand Index evaluated by deriving hard clusters from the
% estimated components
% sortedBasisNum : range of values for which components were estimated

listing = dir(pathDir1);
listing=listing(3:end) ;
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

ARI=zeros(numDifBases,1);

for exp=1:numDifBases
    disp([ num2str(exp) '/' num2str(numDifBases)])
   
    resSplit1 = load([pathDir1 '/NumBases' num2str(sortedBasisNum(exp)) '/OPNMF/ResultsExtractBases.mat']) ;
    resSplit2 = load([pathDir2 '/NumBases' num2str(sortedBasisNum(exp)) '/OPNMF/ResultsExtractBases.mat']) ;
    
    % normalize to unit norm
    wlen1 = sqrt(sum((resSplit1.B).^2)) ;
    wlen2 = sqrt(sum((resSplit2.B).^2)) ;    
    
    if any(wlen1==0)
        wlen1(wlen1==0) = 1;
    end
    W1 = bsxfun(@times,resSplit1.B,1./wlen1) ;
   
    if any(wlen2==0)
        wlen2(wlen2==0) = 1;
    end
    W2 = bsxfun(@times,resSplit2.B,1./wlen2) ;
    
    % calculate inner products
    inner_product = W1'*W2 ;
    
    % take a distance
    dist = 2*(1 - inner_product) ;
    
    % find correspondences
    [Matching,~] = Hungarian(dist);
    [~,idx_hug1]=max(Matching,[],2);
    
    % overlap - hungarian
    overlap{exp} = zeros(length(wlen1),1) ;
    for b=1:length(wlen1)
        overlap{exp}(b) = inner_product(b,idx_hug1(b));
    end
    
    % overlap with best
    overlap_best{exp} = max(inner_product,[],2) ;
    
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
    ARI(exp) = clustering_adjustedRand_fast(clustering1,clustering2);
    
end

meanInner=cellfun(@(x) mean(x),overlap,'UniformOutput', false);
medianInner=cellfun(@(x) median(x),overlap,'UniformOutput', false);

% save results
save([outputDir '/reproducibilityResults.mat'],'meanInner','medianInner','ARI','sortedBasisNum');

%stdInner=cellfun(@(x) std(x),overlap,'UniformOutput', false);
figure;plot(sortedBasisNum,cell2mat(meanInner),'b','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel({'Split-sample reproducibility';'(mean inner product)'},'fontsize',12)
set(gca,'fontsize',12)
saveas(gcf,[outputDir '/MeanInnerProductReproducibility.fig'])
saveas(gcf,[outputDir '/MeanInnerProductReproducibility.png'])

figure;plot(sortedBasisNum,cell2mat(medianInner),'b','LineWidth',2)
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

