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

fid=fopen(list,'r');
if (fid == -1)
    error(['extractBases:loadDataFromList ','Can not open ' list ' file.']);
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
    [status,scratchDir]=system('echo ${SBIA_TMPDIR}');
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
        error(['extractBases:loadDataFromList ','Can not find process id.']);
    end
    % make temp directory to smooth images
    tmpDirName=[scratchDir '/smoothImg_NumBases_' num2str(param.numBases) '_PID_' num2str(pid) ] ;
    
    for tt=1:3
        % if directory exists I will try to create a unique name 3 times
        if(exist(tmpDirName,'dir'))
            tmpDirName = [tmpDirName '_' num2str(floor(rand(1)*10000)) ] ;
        end
    end
    if(~exist(tmpDirName,'dir'))
        status=mkdir(tmpDirName);
        if(status ~= 1)
            error(['extractBases:loadDataFromList ','Can not create temporary directory ' tmpDirName ' to store smoothed images.']);
        end
    else
        error(['extractBases:loadDataFromList ','Someone else is writing in temporary directory ' tmpDirName ]);
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
        error(['extractBases:loadDataFromList ','Mask size does not agree with image size !!']);
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
for i=1:count
    
    % added print statement - Tom
    disp(sprintf("Loading image (%d/%d)...", i, count));
    
    if(~isempty(param.smooth) && param.smooth ~= 0)
        smooth_command = ['3dmerge -1blur_fwhm ' num2str(param.smooth) ' -prefix ' tmpDirName '/smoothed_image_' num2str(i) '.nii.gz ' datafullpath{i} ] ;
        [status,~] = system(smooth_command) ;
        if(status ~= 0)
            error(['extractBases:loadDataFromList ','Can not smooth image ' datafullpath{i} '.']);
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
end
data.nii = nii ;

% % if data were smoothed, remove directory
if(~isempty(param.smooth) && param.smooth ~= 0)
    [status,~,~]=rmdir(tmpDirName);
    if(status ~= 1)
        error(['extractBases:loadDataFromList ','Can not remove folder ' tmpDirName ]);
    end
end
