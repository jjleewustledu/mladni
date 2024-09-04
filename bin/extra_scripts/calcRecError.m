function RecError = calcRecError(X, resultDir, outputDir)

% Function that calculates the reconstruction error given data matrix X and the non-negative matrix factorization results saved in the resultDir directory. The function additionally saves figures that plot the reconstruction error, the gradient of the reconstruction error, and the percentage of improvement as a function of the number of components

% load results and calculate reconstruction error
% We assume that the results directory is organized as follows: there is a folder for every solution (i.e., NumBases${K}, wehere K is the number of components for the solution), and inside the folder there is .mat file (ResultsExtractBases.mat) that contains the matrices W and H that were estimated by the non negative matrix factorization

listing = dir(resultDir);
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

RecError=zeros(numDifBases,1);
for b=1:numDifBases
    disp(b/numDifBases)
    load([resultDir '/NumBases' num2str(sortedBasisNum(b)) '/OPNMF/ResultsExtractBases.mat'])  
    Est = B*C ;
    RecError(b) = norm(X-Est,'fro') ;    
    clear B C
end

% make figures
% 1) reconstruction error
figure;plot(sortedBasisNum,RecError,'r','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel('Reconstruction error','fontsize',12)
xlim([sortedBasisNum(1) sortedBasisNum(end)])
set(gca,'fontsize',12)
saveas(gcf,[outputDir 'RecError.fig'])
saveas(gcf,[outputDir 'RecError.png'])

% 2) gradient of reconstruction error
figure;plot(sortedBasisNum(2:end),diff(RecError),'r','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel('Gradient of reconstruction error','fontsize',12)
xlim([sortedBasisNum(1) sortedBasisNum(end)])
set(gca,'fontsize',12)
saveas(gcf,[outputDir 'gradientRecError.fig'])
saveas(gcf,[outputDir 'gradientRecError.png'])

% 3) Percentage of improvement over range of components used
figure;plot(sortedBasisNum,abs(RecError-RecError(1))./abs(RecError(1)-RecError(end)),'r','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel('Percentage of improvement over range of components used','fontsize',12)
xlim([sortedBasisNum(1) sortedBasisNum(end)])
set(gca,'fontsize',12)
saveas(gcf,[outputDir 'percentageImprovementRecError.fig'])
saveas(gcf,[outputDir 'percentageImprovementRecError.png'])

close all
