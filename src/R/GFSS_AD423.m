% # GFSS AD423 
% 
% GFSS AD423 is a MATLAB script to generate GFSS eigenbrain images from a matrix (included with package) of preprocessed FDG-PET images from 423 subjects across the Alzheimer’s disease spectrum as described in:  
% 
% Jones et al., Patterns of neurodegeneration in dementia reflect a global functional state space 
% medRxiv 2020.11.09.20228676; doi: https://doi.org/10.1101/2020.11.09.20228676 
% 
% ## System requirements 
% 
%     OS: Linux, Windows, macOS 
%     Matlab v9.4 
%     Statistics and Machine Learning Toolbox v11.3 
%     Parallel Computing Toolbox v6.12 
%     MATLAB Distributed Computing Server v6.12 
%     Jimmy Shen (2021). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), MATLAB Central File Exchange. Retrieved April 28, 2021. 
% 
% ## Installation 
% 
% Add the GFSS_AD423 folder to your MATLAB file path (install <1 sec) 
% 
% ## Usage 
% 
% cd to GFSS_AD423 and run the script GFSS_AD423.m 
% 
% Run time < 1 min 
%  
% 
% ## Output 
% 
% Ten NIfTI format eigenbrain images in custom standard space in the folder ./EigenBrains 
% 
% MNI152 space version of these images are available at: https://neurovault.org/collections/AXJZMEAY/ 
% 
% ## License 
% 
% [MIT](https://choosealicense.com/licenses/mit/) 
% 
% David Thomas Jones
% Jones.David@mayo.edu
% 4/28/2021


% set current directory to 'GFSS_AD423' 
EBoutputdir = './EigenBrains/'; % output directory for eigenbrains
load('RASimg.mat') % load subject by voxel matrix for 423 subjects on the AD spectrum preprocessed as described in the manuscript
N = size(RASimg,1); % number of subjects 
load('fpnii.mat') % load header information for standard space 
[i, j, k] = size(fpnii.img);% standard space dimensions 
RASimg_NZ = (RASimg(:,mean(RASimg)>0));% eliminates zero elements in the matrix outside mask
VoxMed = median(RASimg_NZ,1);% median in each voxel for standardization  
VoxIQR = iqr(RASimg_NZ);% IQR in each voxel for standardization 

% non-parametric voxelwsie standardization 
RASimg_NZ_standard = zeros(N,size(RASimg_NZ,2));% preallocate
parfor ind=1:N
    RASimg_NZ_standard(ind,:)=(RASimg_NZ(ind,:)-VoxMed)./VoxIQR;
end

A = RASimg_NZ_standard'; % transpose to voxel by subject matrix
A = A - repmat(mean(A),size(A,1),1); % subject-wise centering 
C = A'*A; % subject by subject covarience 
[U,S,E] = svd(C); % SVD

% scree plot of eigenvalues
EXP = zeros(N,1);% preallocate
parfor ind = 1:N
    EXP(ind) = (S(ind,ind)/trace(S))*100;
end
subplot(1,2,1); plot(EXP(1:end),'o','MarkerFaceColor', 'b'); xlabel('Eigenvalue Rank'); ylabel('Proportion of Variance Explained'); % scree plot
hold on; plot(EXP(1:10),'o','MarkerFaceColor', 'r'); % highlight first 10 in red 
subplot(1,2,2); plot(cumsum(EXP(1:end)),'o','MarkerFaceColor', 'b'); xlabel('Eigenvalue Rank'); ylabel('Total Variance Explained')% plot cumulative variance explained 
hold on; plot(cumsum(EXP(1:10)),'o','MarkerFaceColor', 'r');% highlight first 10 in red

size_E = 10; % 10 is 95% and reasonable scree
Ee = E(:,1:1:size_E); % eigenvalues for GFSS
EB = A * Ee; % eigenbrains for GFSS

% Convert EBs back to image space and save out images as nifti files 
d = zeros(i*j*k,size_E);
d((mean(RASimg)>0),:) = EB(:,:);
EB4D = zeros(size_E, i, j, k);
for ind = 1:size_E
    eb = reshape(d(:,ind), i, j, k);
    fpnii.img = eb; 
    save_nii(fpnii,[EBoutputdir 'AD423_EB_' num2str(ind,'%03.f') '.nii'])
    EB4D(ind,:,:,:) = eb;
end

% Calculate eigenvalues from preprocessed images using GFSS EB images
ev = zeros(N, size_E);
for ind1 = 1:N  
  for ind2 = 1:size_E 
    ev(ind1,ind2) =   sum(A(:,ind1).* EB(:,ind2))./S(ind2,ind2);
  end
end
