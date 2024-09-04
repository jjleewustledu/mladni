clear
clc
close all

% homePath = '/home/sotirasa/sbia_server/home/' ;
homePath = '/sbia/home/sotirasa/' ;

% populating path
% addpath(genpath('../../'));

mfMeth = {'OPNMF','PCA','Jade'} ;

% load data -> we need to estimate the age
param.isList = 1 ;
param.downSample = 1 ;
dataInput = [ homePath 'Data/BLSA/Lists/BLSA_baseline_pre2010_smoothed_GM.lst' ] ;
data = loadData(dataInput,param,[]) ;
meanX = mean(data.X,2);
X = data.X(meanX>0,:); clear data

% hard coded results path
resultsPath = [ homePath 'comp_space/Experiments/Final/BLSA' ];
listing = dir(resultsPath);
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
        dataPath=[resultsPath '/' listing(idx(i)).name '/' mfMeth{j} '/ResultsExtractBases.mat'];
        % loading B (bases) and C (loading coefficients) - we do not
        % need B
        load(dataPath);
        if(strcmp(mfMeth{j},'PCA'))
            Est = B(meanX>0,:)*C + meanX(meanX>0)*ones(1,size(C,2)); clear B C
        else
            Est = B(meanX>0,:)*C ; clear B C
        end
        RecError(i,j) = norm(X-Est,'fro') ;
        clear Est
    end
end

disp('End')
