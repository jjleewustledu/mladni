 classdef NMF < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% baseline_cn ~ 22 bases
    %  
    %  Created 23-Jun-2022 13:02:04 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
        
    properties (Constant)
        MAX_NUM_BASES = 40
        N_PATTERNS = 24
    end

    methods (Static)
        function create_montage(varargin)
            %  Args:
            %      path (folder): e.g., .,
            %                          /path/to/nmf_dataset, 
            %                          /path/to/nmf_dataset/NumBases20, 
            %                          /path/to/nmf_dataset/NumBases20/OPNMF/niiImg, 
            
            import mladni.NMF;

            ip = inputParser;
            addOptional(ip, 'path', pwd, @isfolder);
            addParameter(ip, 'noclobber', false, @islogical)
            parse(ip, varargin{:});
            ipr = ip.Results;
            if strcmpi(mybasename(ipr.path), 'niiImg') % /path/to/niiImg
                pwd0 = pushd(ipr.path);
                atl = mlfourd.ImagingContext2( ...
                    fullfile(getenv('FSLDIR'), 'data', 'standard', 'MNI152_T1_1mm.nii.gz'));
                for g = glob('Basis_*.nii')'
                    try
                        ic = mlfourd.ImagingContext2(g{1});
                        product = strrep(strcat(ic.fqfp, '.png'), 'niiImg', 'Figures');
                        if ~isfile(product) || ...
                                (isfile(product) && ~ipr.noclobber)
                            atl.save_qc(ic);
                            png = strrep(g{1}, '.nii', '.png');
                            movefile(png, '../Figures', 'f');
                        end
                    catch ME
                        fprintf("%s:  %s\n", stackstr(), ME.message)
                    end
                end
                ensuredir('../Figures');
                popd(pwd0);
                return
            end
            if contains(basename(ipr.path), 'NumBases') % /path/to/NumBases*
                pth_to_niiImg = fullfile(ipr.path, 'OPNMF', 'niiImg', '');
                if isfolder(pth_to_niiImg)
                    NMF.create_montage(pth_to_niiImg, 'noclobber', ipr.noclobber);
                end
                return
            end
            if ~isempty(glob(fullfile(ipr.path, 'NumBases*'))) % nmf_dataset
                for g = glob(fullfile(ipr.path, 'NumBases*'))'
                    NMF.create_montage(g{1}, 'noclobber', ipr.noclobber);
                end
                return
            end
            warning('mladni:NotImplementedWarning', 'NMF.create_montage: nothing to be done');
        end
        function create_mean_imaging(varargin)
            %% creates mean, var, std
            %  Args:
            %       csv {mustBeFile}
            %       outputDir {mustBeFolder}

            ip = inputParser;
            addRequired(ip, 'csv', @isfile);
            addParameter(ip, 'outputDir', pwd, @isfolder);
            parse(ip, varargin{:});
            ipr = ip.Results;            
            
            csv = readlines(ipr.csv); % string vector

            % create mean
            ic = mlfourd.ImagingContext2(char(csv(1)));
            ic.selectMatlabTool();
            re = regexp(ic.fileprefix, '(?<sub>sub-\d{3}S\d{4}_)(?<ses>ses-\d+_)(?<descrip>\S+)', 'names');
            fp = sprintf('all_%s_mean', re.descrip);
            for i = 2:length(csv)
                try
                    ic = ic + mlfourd.ImagingContext2(char(csv(i)));
                    ic.fileprefix = fp;
                catch
                end
            end
            ic = ic/length(csv);
            ic.fileprefix = fp;
            ic.filepath = ipr.outputDir;
            ic.save();
            
            % create variance
            ic1 = mlfourd.ImagingContext2(char(csv(1)));
            ic1.selectMatlabTool();
            ic1 = (ic1 - ic).^2;
            re1 = regexp(ic.fileprefix, '(?<descrip>\S+)_mean', 'names');
            fp1 = sprintf('%s_variance', re1.descrip);
            for i = 2:length(csv)
                try
                    ic1 = ic1 + (mlfourd.ImagingContext2(char(csv(i))) - ic).^2;
                    ic1.fileprefix = fp1;
                catch
                end
            end
            ic1 = ic1/(length(csv) - 1);
            ic1.fileprefix = fp1;
            ic1.filepath = ipr.outputDir;
            ic1.save();

            % create std
            ic2 = copy(ic1);
            ic2 = sqrt(ic2);
            ic2.fileprefix = strrep(ic1.fileprefix, "_variance", "_std");
            ic2.save();
        end
        function create_median_imaging(varargin)
            %% creates median, iqr
            %  Args:
            %       csv {mustBeFile}
            %       outputDir {mustBeFolder}

            ip = inputParser;
            addRequired(ip, 'csv', @isfile);
            addParameter(ip, 'outputDir', pwd, @isfolder);
            parse(ip, varargin{:});
            ipr = ip.Results;                        

            csv = readlines(ipr.csv); % string vector

            % create median
            ifc = mlfourd.ImagingFormatContext2(char(csv(1)));
            re = regexp(ifc.fileprefix, '(?<sub>sub-\d{3}S\d{4}_)(?<ses>ses-\d+_)(?<descrip>\S+)', 'names');
            fp = sprintf('all_%s_median', re.descrip);
            size3d = size(ifc.img);
            for i = 2:length(csv)
                try
                    ifc1 = mlfourd.ImagingFormatContext2(char(csv(i)));
                    if isempty(ifc1.img)
                        ifc.img(:,:,:,i) = nan(size3d);
                    else
                        ifc.img(:,:,:,i) = ifc1.img;
                    end
                catch
                    ifc.img(:,:,:,i) = nan(size3d);
                end
            end
            ifc.img = median(ifc.img, 4, 'omitnan');
            ifc.fileprefix = fp;
            ifc.filepath = ipr.outputDir;
            ifc.save();     

            % create iqr
            ifc = mlfourd.ImagingFormatContext2(char(csv(1)));
            re = regexp(ifc.fileprefix, '(?<sub>sub-\d{3}S\d{4}_)(?<ses>ses-\d+_)(?<descrip>\S+)', 'names');
            fp = sprintf('all_%s_iqr', re.descrip);
            size3d = size(ifc.img);
            for i = 2:length(csv)
                try
                    ifc1 = mlfourd.ImagingFormatContext2(char(csv(i)));
                    if isempty(ifc1.img)
                        ifc.img(:,:,:,i) = nan(size3d);
                    else
                        ifc.img(:,:,:,i) = ifc1.img;
                    end
                catch
                    ifc.img(:,:,:,i) = nan(size3d);
                end
            end
            ifc.img = iqr(ifc.img, 4); % safe for nans
            ifc.fileprefix = fp;
            ifc.filepath = ipr.outputDir;
            ifc.save();   
        end        
        function csv_out = filter_missing_images(varargin)
            %% reads csv, ignores missing rows, adjusts SINGULARITY_HOME as needed, writes csv_out
            
            ip = inputParser;
            addRequired(ip, 'csv', @isfile);
            addParameter(ip, 'csv_out', '', @istext)
            addParameter(ip, 'singularity_home', getenv('SINGULARITY_HOME'), @isfolder);
            parse(ip, varargin{:});
            ipr = ip.Results;  
            if isempty(ipr.csv_out)
                [pth,fp,x] = fileparts(ipr.csv);
                ipr.csv_out = fullfile(pth, strcat(fp, '_', x));
            end

            %fprintf('mladni.NMF.filter_missing_images():\n')
            
            lines = readlines(ipr.csv); % string vector
            lines1 = {};
            Nmissing = 0;
            for li = 1:length(lines)
                line_ = lines{li};
                if isempty(line_)
                    break
                end
                if ~contains(line_, ipr.singularity_home)
                    ss = strsplit(line_, 'Singularity');
                    line_ = strcat(ipr.singularity_home, ss{2});
                end
                if isfile(line_)
                    lines1 = [lines1; line_];
                else
                    fprintf('\tmissing %s\n', line_);
                    Nmissing = Nmissing + 1;
                end
            end

            writetable(table(lines1), ipr.csv_out, 'WriteVariableNames', false);
            if Nmissing > 0
                fprintf('\tdiscovered %i missing files\n', Nmissing);
            end
            fprintf('\tread %s\n', ipr.csv);
            fprintf('\twrote %s\n', ipr.csv_out);
            csv_out = ipr.csv_out;
        end      
        function diagnose_reproducibility()
            this = mladni.NMF(nmfDataset='baseline_cn');
            cd(this.outputDir);            
            this.evaluateRepeatedReproducibility2('cn')
        end                
        function split_nifti_files(csv)
            %% DEPRECATED; see mladni.Adni.anticlust_cn_repeat

            [pth,fp,x] = myfileparts(csv);
            if isempty(pth)
                pth = pwd;
            end
            assert(isfolder(pth));
            assert(strcmp(x, '.csv'));
            pth1 = strcat(pth, '_split1');
            ensuredir(pth1)
            pth2 = strcat(pth, '_split2');
            ensuredir(pth2)

            % read table to split
            tbl = readtable(csv, 'ReadVariableNames', false, 'Delimiter', ' ');

            % define splits
            s1 = rand(size(tbl.Var1));
            s1 = s1 > 0.5;
            s2 = ~s1;

            % make splits
            Var1 = tbl.Var1(s1);
            csv1 = fullfile(pth1, strcat(fp, x));
            writetable(table(Var1), csv1, 'WriteVariableNames', false);

            Var2 = tbl.Var1(s2);
            csv2 = fullfile(pth2, strcat(fp, x));
            writetable(table(Var2), csv2, 'WriteVariableNames', false);
        end
        
        %% CHPC functions

        function propcluster()
            %% PROPCLUSTER.
            %  Currently, the script asks for 24Gb RAM, 4 CPUs and 24h of compute time. Asking for more CPUs will make 
            %  the executable finish faster. The more resources one asks, the more time it will probably take for the 
            %  job to start. If an executable needs more resources than the ones requested, it will fail. I am not sure 
            %  about the exact requirements here as this depends both on the number of samples and the dimensionality of 
            %  data. I am almost certain that the nmf computations will take more than 24h (especially for the higher 
            %  number of components). However, I would be sKeptical about requesting more time as it seems that the 
            %  queueing system heavily penalizes longer jobs (small number of max running jobs). Note that the code 
            %  saves intermediate results and can restart from these intermediate points by simply running the same 
            %  command and pointing to the same output directory.

            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '128000'; 
            c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = '';
            c.AdditionalProperties.WallTime = '24:00:00'; % 24 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function parcall()
            %% DEPRECATED legacy method
            groups = mladni.NMF.groups;
            parfor gi = 1:length(groups)
                try
                    this = mladni.NMF( ...
                        nmfDataset=sprintf('baseline_%s', groups{gi}));
                    cd(this.outputDir);
                    call(this);
                    this.evaluateRepeatedReproducibility(groups{gi})
                catch ME
                    handwarning(ME)
                end
            end
        end
        function parcall2()
            %% DEPRECATED legacy method
            targetDatasetDir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG', 'baseline_cn');
            assert(isfolder(targetDatasetDir))

            b = 2:2:mladni.NMF.MAX_NUM_BASES;
            parfor bi = 1:length(b)
                groups = mladni.NMF.groups;
                for gi = 1:length(groups)
                    try
                        this = mladni.NMF( ...
                            selectedNumBases=b(bi), ...
                            nmfDataset=sprintf('baseline_%s', groups{gi}));
                        cd(this.outputDir);
                        call2(this, targetDatasetDir, 'on_cn');
                    catch ME
                        handwarning(ME)
                    end
                end
            end

            %% gather component files on local machine
            % for b in baseline_*; do
            % pushd $b
            % for n in NumBases*; do
            % pushd $n
            % rsync -ra login3.chpc.wustl.edu:/scratch/jjlee/Singularity/ADNI/NMF_FDG/$b/$n/components .
            % popd; done
            % popd; done
        end        
        function par_create_montage(gid)
            arguments
                gid {mustBeScalarOrEmpty} = 1
            end

            if gid
                groups = mladni.NMF.groups;
                this = mladni.NMF(nmfDataset=sprintf('baseline_%s', groups{gid}));
                cd(this.outputDir);
                this.create_montage(this.outputDir, noclobber=true);
            else
                groups = mladni.NMF.groups;
                parfor gi = 1:length(groups)
                    this = mladni.NMF(nmfDataset=sprintf('baseline_%s', groups{gi}));
                    cd(this.outputDir);
                    this.create_montage(this.outputDir, noclobber=true);
                end
            end
        end
        function getDebugLog(j,c)
            try
                c.getDebugLog(j)
            catch
                c.getDebugLog(j.Tasks(end))
            end
        end

        %% Aris' inference, called from other method functions
        
        function calculateComponentWeightedAverageNIFTI(dataList,resultsDir,numBases,outPath)
            %% loops through NumBases*, then accesses all /OPNMF/niiImg/Basis_*.nii to project files in dataList to Basis_*.nii
            % dataList: .csv file containing the images for which the coefficients need
            %           to be calculated. Full path for every image file is given in
            %           every line
            % resultsDir: directory where the NMF results are (i.e., the level where
            %             the NumBases folder is placed)
            % numBases: determines the solution for which one wants to calculate
            %           subject coefficients
            % outPath: fq-filenamne for saving output
            
            import mlniftitools.*;
            resultsDir = convertStringsToChars(resultsDir);
            outPath = convertStringsToChars(outPath);
            
            for b=1:length(numBases)
                disp(numBases(b))
                dataPath=[resultsDir '/NumBases' num2str(numBases(b)) '/OPNMF/niiImg/'];
                
                % loading estimated non-negative components
                listing = dir(dataPath);
                hh =cellfun(@(x) ( (strfind(x,'Basis')==1)  ),{listing(:).name},'UniformOutput',false) ;
                listing=listing(cellfun(@(x) (~isempty(x)&&(x==1)),hh));
                
                if(length(listing)~=numBases(b))
                    error(['I cannot find ' num2str(numBases(b)) ' basis images in the folder ' dataPath ])
                end
                
                for i=1:numBases(b)
                    nii = load_untouch_nii([dataPath listing(i).name]);
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
                
                % read list
                fid=fopen(dataList,'r');
                if (fid == -1)
                    error(['extractBases:calculateComponentWeightedAverage ','Can not open ' list ' file.']);
                end
                datafullpath = textscan(fid,'%s\n');
                fclose(fid);
                
                datafullpath = datafullpath{1,1} ;
                datafullpath = cellstr(datafullpath) ;
                count = numel(datafullpath);
                
                fid = fopen(outPath,'w');
                frmtWrite='%s,';
                frmtWrite=[frmtWrite repmat('%f,',1,numBases(b)-1)]; %#ok<*AGROW>
                frmtWrite=[frmtWrite '%f\n'];
                
                wA = zeros(count,numBases) ;
                  
                %print statement added by tom
                disp(sprintf('calculateComponentWeightedAverageNIFTI: loading %d files', count)); %#ok<DSPS>
                for i=1:count                    
                    nii = load_untouch_nii(datafullpath{i});
                    wA(i,:) = double(nii.img(:)')*nB;
                    fprintf(fid,frmtWrite,datafullpath{i},wA(i,:)');
                end
                fclose(fid);
            end
        end
        function calculateSelectedComponentWeightedAverageNIFTI(dataLists,resultsDir,numBases,outFile)
            %% selects single folder NumBases\d, then accesses all /OPNMF/niiImg/Basis_*.nii to project files in dataList to Basis_*.nii
            % dataList: .csv file containing the images for which the coefficients need
            %           to be calculated. Full path for every image file is given in
            %           every line
            % resultsDir: directory where the NMF results are (i.e., the level where
            %             the NumBases folder is placed)
            % numBases: determines the solution for which one wants to calculate
            %           subject coefficients
            % outFile: CSV-file saving output
            
            import mlniftitools.*;
            resultsDir = convertStringsToChars(resultsDir);
            outFile = convertStringsToChars(outFile);
        
            % works with only single folder NumBases\d found from reproducibility analysis

            fprintf("%s: numBases->%g\n", stackstr(), numBases)
            dataPath=[resultsDir '/NumBases' num2str(numBases) '/OPNMF/niiImg/'];
            
            % loading estimated non-negative components
            listing = dir(dataPath);
            hh =cellfun(@(x) ( (strfind(x,'Basis')==1)  ),{listing(:).name},'UniformOutput',false) ;
            listing=listing(cellfun(@(x) (~isempty(x)&&(x==1)),hh));
            
            if(length(listing)~=numBases)
                error(['I cannot find ' num2str(numBases) ' basis images in the folder ' dataPath ])
            end
            
            for i=1:numBases
                nii = load_untouch_nii([dataPath listing(i).name]);
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
        function [RecError,X] = calcBLSARecError(inFiles, resultDir)

            resultDir = convertStringsToChars(resultDir);

            % populating path
            % addpath(genpath('../../'));

            mfMeth = {'OPNMF','PCA','Jade'} ;

            % load data -> we need to estimate the age
            param.isList = 1 ;
            param.downSample = 1 ;
            data = loadData(inFiles,param,[]) ;
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
        function RecError = calcRecError(X, resultDir, resultsDir, selectVoxels)
            
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
                load([resultDir '/NumBases' num2str(sortedBasisNum(b)) '/OPNMF/ResultsExtractBases.mat']) %#ok<LOAD>
                Est = B*C ;
                assert(size(Est,1) == size(selectVoxels,1), stackstr()) %% JJL's BUG CHECK
                Est = Est(selectVoxels,:); %% JJL
                RecError(b) = norm(X-Est,'fro') ;
                clear B C
            end
            
            % make figures
            % 1) reconstruction error
            mladni.NMF.plotForPub(sortedBasisNum, RecError, ...
                xlab="Number of patterns", ...
                ylab="Reconstruction error", ...
                fileprefix=fullfile(resultsDir, "RecError"));
            
            % 2) gradient of reconstruction error
            mladni.NMF.plotForPub(sortedBasisNum(2:end), diff(RecError), ...
                xlab="Number of patterns", ...
                ylab="Gradient of reconstruction error", ...
                fileprefix=fullfile(resultsDir, "gradientRecError"));
            
            % 3) Percentage of improvement over range of components used
            improvement = 100*abs(RecError-RecError(1))./abs(RecError(1)-RecError(end));
            mladni.NMF.plotForPub(sortedBasisNum, improvement, ...
                xlab="Number of patterns", ...
                ylab="% improvement with cumulative patterns", ...
                fileprefix=fullfile(resultsDir, "percentageImprovementRecError"));
            
            %close all
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
        function tf = components_are_available()
            tf = false;
        end        
        function evaluateRepeatedReproducibility(subgroup, N, K, saveGcf)
            arguments
                subgroup {mustBeTextScalar} = 'cn'
                N double {mustBeInteger} = 50 % # repetitions
                K double {mustBeInteger} = 20 % # bases checked by VolBin/*
                saveGcf logical = true
            end

            home = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
            subgroups = { ...
                'cn', 'preclinical', ...
                'cdr_0p5_apos', ...
                'cdr_gt_0_aneg', 'cdr_gt_0p5_apos'};
            assert(any(strcmp(subgroup, subgroups)));
            
            outputDir_ = fullfile(home, sprintf('baseline_%s', subgroup));
            ensuredir(outputDir_);
            ARI = zeros(K, N);
            overlap = cell(K, N);
            for rep = 1:N
                try
                    pathDir1_ = fullfile(home, sprintf('baseline_%s_repA%i', subgroup, rep));
                    pathDir2_ = strrep(pathDir1_, 'repA', 'repB');
                    [~,~,A_,o_] = mladni.NMF.evaluateReproducibility( ...
                        pathDir1_, pathDir2_, outputDir_, false);
                    ARI(:,rep) = A_;
                    for ik = 1:K
                        try
                            overlap{ik,rep} = o_{ik};
                        catch
                        end
                    end
                catch ME
                    handwarning(ME)
                end
            end

            mu_ARI = mean(ARI, 2)';
            sigma_ARI = std(ARI, 0, 2)';

            mu_overlap = cell(1, K);
            sigma_overlap = cell(1, K);
            for basis = 1:K
                try
                    accum = [overlap{basis, :}];
                    mu_overlap{basis} = mean(accum, "all");
                    sigma_overlap{basis} = std(accum, 0, "all");
                catch 
                end
            end
            mu_overlap = cell2mat(mu_overlap);
            sigma_overlap = cell2mat(sigma_overlap);
            save(fullfile(outputDir_, 'evaluateRepeatedReproducibility.mat'), ...
                'ARI', 'overlap', 'mu_ARI', 'sigma_ARI', 'mu_overlap', 'sigma_overlap');

            if saveGcf
                sortedBasisNum = 2:2:2*K;
                xconf = [sortedBasisNum sortedBasisNum(end:-1:1)];

                figure
                yconf = [mu_overlap+sigma_overlap mu_overlap(end:-1:1)-sigma_overlap];
                p = fill(xconf,yconf, 'red');
                p.FaceColor = [0.6 0.8 1];
                p.EdgeColor = 'none';
                hold on
                plot(sortedBasisNum, mu_overlap, 'bo-', 'LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',[0 0 0])
                hold off
                xlabel('Number of components','fontsize',12)
                ylabel({'Split-sample reproducibility';'mean +/- std of inner product'},'fontsize',12)
                title(strcat('baseline_', subgroup), Interpreter="none")
                set(gca,'fontsize',12)
                saveas(gcf, fullfile(outputDir_, 'innerProductReproducibility.fig'))
                saveas(gcf, fullfile(outputDir_, 'innerProductReproducibility.png'))
                
                figure
                yconf = [mu_ARI+sigma_ARI mu_ARI(end:-1:1)-sigma_ARI];
                p = fill(xconf,yconf, 'red');
                p.FaceColor = [0.6 0.8 1];
                p.EdgeColor = 'none';                
                hold on
                plot(sortedBasisNum, mu_ARI, 'bo-', 'LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',[0 0 0])
                hold off
                xlabel('Number of components','fontsize',12)
                ylabel({'Split-sample reproducibility' ;'mean +/- std of Adjusted Rand Index'},'fontsize',12)
                title(strcat('baseline_', subgroup), Interpreter="none")
                set(gca,'fontsize',12)
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility.fig'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility.png'))
            end
        end
        function evaluateRepeatedReproducibility2(subgroup, N, K, saveGcf, opts)
            arguments
                subgroup {mustBeTextScalar} = 'cn'
                N double {mustBeInteger} = 50 % # repetitions
                K double {mustBeInteger} = 20 % # bases checked by VolBin/*
                saveGcf logical = true
                opts.fracx double = 0.5
                opts.Npx double = 3400
                opts.fracy double = 0.5
                opts.Npy double = 1440
            end

            import mladni.Jones2022.al_goodplot

            home = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
            subgroups = { ...
                'cn', 'preclinical', ...
                'cdr_0p5_apos', ...
                'cdr_gt_0_aneg', 'cdr_gt_0p5_apos'};
            assert(any(strcmp(subgroup, subgroups)));
            
            outputDir_ = fullfile(home, sprintf('baseline_%s', subgroup), 'results');
            ensuredir(outputDir_);
            ARI = nan(K, N);
            overlap = cell(K, N);

            if ~isfile(fullfile(outputDir_, 'evaluateRepeatedReproducibility2.mat'))
                parfor rep = 1:N
                    try
                        pathDir1_ = fullfile(home, sprintf('baseline_%s_repA%i', subgroup, rep));
                        pathDir2_ = strrep(pathDir1_, 'repA', 'repB');
                        [~,~,A_,o_] = mladni.NMF.evaluateReproducibility( ...
                            pathDir1_, pathDir2_, outputDir_, false);
                        ARI(:,rep) = A_;
                        for ik = 1:K
                            try
                                overlap{ik,rep} = o_{ik};
                            catch %#ok<CTCH>
                            end
                        end
                    catch ME
                        handwarning(ME)
                    end
                end
                ARI = mladni.NMF.removeColsWithNans(ARI);
                overlap = mladni.NMF.removeColsWithNans(overlap);
                save(fullfile(outputDir_, 'evaluateRepeatedReproducibility2.mat'), ...
                    'ARI', 'overlap');
            else
                ld = load(fullfile(outputDir_, 'evaluateRepeatedReproducibility2.mat'));
                ARI = ld.ARI;
                overlap = ld.overlap;
            end

            if saveGcf
                sortedBasisNum = 2:2:2*K;
                sortedBasisNames = cellfun(@num2str, num2cell(2:2:2*K), UniformOutput=false);

                meanInner=cell2mat(cellfun(@(x) mean(x),overlap,'UniformOutput', false));
                medianInner=cell2mat(cellfun(@(x) median(x),overlap,'UniformOutput', false));           

                figure
                al_goodplot(sortedBasisNum, sortedBasisNames, ARI');
                xlabel('Number of patterns','fontsize',20)
                ylabel({'Split-sample reproducibility' ;'(Adjusted Rand Index)'},'fontsize',20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat2.fig'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat2.png'))
                saveas(gcf, fullfile(outputDir_, 'AdjustedRandIndexReproducibility_repeat2.svg'))

                figure
                al_goodplot(sortedBasisNum, sortedBasisNames, meanInner');
                xlabel('Number of patterns','fontsize',20)
                ylabel({'Split-sample reproducibility';'(mean inner product)'},'fontsize',20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat2.fig'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat2.png'))
                saveas(gcf, fullfile(outputDir_, 'MeanInnerProductReproducibility_repeat2.svg'))

                figure
                al_goodplot(sortedBasisNum, sortedBasisNames, medianInner');
                xlabel('Number of patterns','fontsize',20)
                ylabel({'Split-sample reproducibility';'(median inner product)'},'fontsize',20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat2.fig'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat2.png'))
                saveas(gcf, fullfile(outputDir_, 'MedianInnerProductReproducibility_repeat2.svg'))

                ARI_snr = median(ARI', 1)./iqr(ARI', 1); %#ok<UDIM>
                mladni.NMF.plotForPub(sortedBasisNum, ARI_snr, ...
                    xlab="Number of patterns", ...
                    ylab=["Split-sample reproducibility"; "(adjusted Rand index, median/iqr)"], ...
                    fileprefix=fullfile(outputDir_, "MedianIqrARIReproducibility_repeat2"));
            end
        end
        function [meanInner,medianInner,ARI,overlap,sortedBasisNum] = evaluateReproducibility( ...
                pathDir1, pathDir2, outputDir, saveGcf)
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
            %
            % outputs
            % meanInner : mean value of the inner product between matched components
            % medianInner : median value of the inner product between matched components
            % ARI : adjusted Rand Index evaluated by deriving hard clusters from the estimated components
            % overlap : double, size ~ size ~ {1,numDifBases}[length(wlen1) in 2:2:40, 1]
            % sortedBasisNum : size ~ [numDifBases,1]
            
            arguments
                pathDir1 char {mustBeTextScalar,mustBeFolder}
                pathDir2 char {mustBeTextScalar,mustBeFolder}
                outputDir char {mustBeTextScalar}
                saveGcf logical = true
            end
            pathDir1 = convertStringsToChars(pathDir1);
            pathDir2 = convertStringsToChars(pathDir2);
            outputDir = convertStringsToChars(outputDir);

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
            
            for exp_=1:numDifBases
                disp([ num2str(exp_) '/' num2str(numDifBases)])
               
                resSplit1 = load([pathDir1 '/NumBases' num2str(sortedBasisNum(exp_)) '/OPNMF/ResultsExtractBases.mat']) ;
                resSplit2 = load([pathDir2 '/NumBases' num2str(sortedBasisNum(exp_)) '/OPNMF/ResultsExtractBases.mat']) ;
                
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
                
            end
            
            meanInner=cell2mat(cellfun(@(x) mean(x),overlap,'UniformOutput', false));
            medianInner=cell2mat(cellfun(@(x) median(x),overlap,'UniformOutput', false));
            
            % save results
            save([outputDir '/reproducibilityResults.mat'],'meanInner','medianInner','ARI','sortedBasisNum');
            if saveGcf
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
                    error(['extractBases:loadDataFromList ','Can not find process id.']);
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
            % added print statement - Tom
            disp(sprintf("Loading %d images from %s ...", count, param.componentDir)); %#ok<DSPS>
            for i=1:count
                                
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
                    error(['*:loadDataFromList ','Can not remove folder ' tmpDirName ]);
                end
            end
        end
        function data = loadData(dataInput,param,subsetIdx)
            
            import mlniftitools.*;
            import mladni.NMF.loadDataFromList;
            
            if(param.isList==1)
                disp(['List ' dataInput ]);
                data = loadDataFromList(dataInput,param,subsetIdx) ; % returned data is a struct !
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
        function h = plotForPub(x, y, opts)
            %% PLOTFORPUB provides a consistent plotting style for publications, also saving {.fig,.png,.svg}.
            %  For plot h, h.Position := [1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy]
            %  Args:
            %  x double
            %  y double
            %  opts.xlab {mustBeText} = ""
            %  opts.ylab {mustBeText} = ""
            %  opts.fileprefix {mustBeTextScalar} = stackstr(3)
            % opts.fracx double = 0.5
            % opts.Npx double = 3400
            % opts.fracy double = 0.5
            % opts.Npy double = 1440

            arguments
                x double
                y double
                opts.xlab {mustBeText} = ""
                opts.ylab {mustBeText} = ""
                opts.fileprefix {mustBeTextScalar} = ""
                opts.fracx double = 0.5
                opts.Npx double = 3400
                opts.fracy double = 0.5
                opts.Npy double = 1440
            end

            h = figure;
            plot(ascol(x), ascol(y), ...
                ':_', LineWidth=3, Color=[0.73,0.83,0.96], ...
                MarkerSize=25, MarkerEdgeColor=[0,0,0], MarkerFaceColor=[0,0,0]);
            %xlim([x(1), x(end)]);
            xlabel(opts.xlab,'fontsize',20);
            ylabel(opts.ylab,'fontsize',20);
            set(gca,'fontsize',20);
            set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy]);
            if ~isemptytext(opts.fileprefix)
                saveas(gcf, opts.fileprefix+".fig");
                saveas(gcf, opts.fileprefix+".png");
                saveas(gcf, opts.fileprefix+".svg");
            end
        end
        function X = removeColsWithNans(X)
            if iscell(X)
                X = X(:, all(~isnan(cell2mat(X))));
                return
            end                        
            % remove reps (columns) that have any nan elements
            X = X(:, all(~isnan(X)));
        end
    end
    
    properties
        cache_files
        downSample
        isList
        mask
        mcrroot
        memInGB
        nmfDataset
        numBases
        permute
        repetitions
        selectedNumBases
        smooth
        study_design
        volbin
    end

    properties (Constant)
        groups = { ...
            'cn', 'preclinical', ...
            'cdr_0p5_apos', ...
            'cdr_gt_0p5_apos', ...
            'cdr_gt_0_aneg'}
    end

    properties (Dependent)
        adniDemographics
        componentDir
        groupPrefix
        inFiles
        inFiles2
        isCrossSectional
        outputDir
        targetDatasetDir % baseline_cn that provides all NMF patterns
        Xmat
        X2mat
    end

    methods % GET
        function g = get.adniDemographics(this)
            if ~isempty(this.adniDemographics_)
                g = this.adniDemographics_;
                return
            end
            this.adniDemographics_ = mladni.AdniDemographics(study_design=this.study_design);
            g = this.adniDemographics_;
        end
        function g = get.componentDir(this)
            g = fullfile(this.outputDir, sprintf('NumBases%i', this.selectedNumBases), 'components');
            ensuredir(g);
        end
        function g = get.groupPrefix(this)
            switch convertCharsToStrings(this.study_design)
                case "cross-sectional"
                    g = "baseline_";
                case "longitudinal"
                    g = "all_";
                otherwise
                    error("mladni:ValueError", "%s: this.study_design->%s", stackstr(), this.study_design)
            end            
        end
        function g = get.inFiles(this)
            g = fullfile(this.outputDir, 'nifti_files_mounted.csv');
            if ~isfile(g)
                g0 = strrep(g, '_mounted.csv', '.csv');
                t = readtable(g0, Format='%s', Delimiter=' ', ReadVariableNames=false);
                t.Var1 = strrep(t.Var1, '/scratch/jjlee', '/home/usr/jjlee');
                writetable(t, g, WriteVariableNames=false);
            end
        end
        function g = get.inFiles2(this)
            nmfc = mladni.NMFCovariates();
            g = nmfc.inFiles;
        end
        function g = get.isCrossSectional(this)
            g = contains(this.study_design, 'cross', IgnoreCase=true) && ...
                contains(this.study_design, 'sectional', IgnoreCase=true);
        end
        function g = get.outputDir(this)
            g = fullfile(getenv('ADNI_HOME'), 'NMF_FDG', this.nmfDataset, '');
        end
        function g = get.targetDatasetDir(this)
            g = fullfile(getenv('ADNI_HOME'), 'NMF_FDG', 'baseline_cn');
        end
        function g = get.Xmat(this)
            numBasesFolder = sprintf('NumBases%i', this.selectedNumBases);
            g = fullfile(this.outputDir, numBasesFolder, 'X.mat');
        end
        function g = get.X2mat(this)
            numBasesFolder = sprintf('NumBases%i', this.selectedNumBases);
            g = fullfile(this.outputDir, numBasesFolder, 'X2.mat');  
        end
    end
    
    methods
        function build_for_brainsmash(this)
            %% Prepares intermediates needed for brainsmash inference described in 
            %  test_neurodegeneration2.py/TestNeurodegeneration2.test_build_stats27.
            %  Requires availability of AFNI.3dmaskdump.

            workdir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG', 'baseline_cn', ...
                sprintf('NumBases%i', this.selectedNumBases), 'OPNMF', 'niiImg');
            pwd1 = pushd(workdir);
            parfor idx = 1:this.selectedNumBases
                system(sprintf( ...
                    '3dmaskdump -mask %s -o Basis_%i.txt -noijk Basis_%i.nii', ...
                    fullfile(getenv('ADNI_HOME'), 'VolBin', 'mask.nii.gz'), ...
                    idx, idx));
            end
            popd(pwd1)
        end
        function build_stats_imaging(this, opts)
            %% This apex function calls create_{montage, mean_imaging, median_imaging}.

            arguments
                this mladni.NMF
                opts.inputDir {mustBeFolder} = pwd     % e.g., NMF_FDG/baseline_cn, NMF_FDG/longitudinal_cdr_0.5_apos, ...
                                                       % has nifti_files_mounted.csv visible to glob()
                opts.outputDir {mustBeTextScalar} = "" 
            end
            assert(contains(opts.inputDir, "NMF"))
            if isemptytext(opts.outputDir)
                opts.outputDir = opts.inputDir;
            end
            import mladni.NMF.*
            
            % montage
            if contains(opts.inputDir, "baseline")
                globbed = glob(fullfile(opts.inputDir, sprintf("NumBases%i", this.selectedNumBases)));
                if isempty(globbed)
                    globbed = glob(fullfile(opts.inputDir, "*", sprintf("NumBases%i", this.selectedNumBases)));
                end
                assert(~isempty(globbed))
                this.create_montage(globbed{1});
            end

            % mean, var, std, median, iqr
            globbed = glob(fullfile(opts.inputDir, "nifti_files_mounted.csv"));
            if isempty(globbed)
                globbed = glob(fullfile(opts.inputDir, "*", "nifti_files_mounted.csv"));
            end
            assert(~isempty(globbed))
            this.create_mean_imaging(globbed{1}, outputDir=opts.outputDir);
            this.create_median_imaging(globbed{1}, outputDir=opts.outputDir);            
        end
        function build_surfaces(this)
            pwd0 = pushd(fullfile(this.outputDir, sprintf('NumBases%i', this.selectedNumBases), 'OPNMF', 'niiImg'));
            mlshimony.Workbench.vol2surf('Basis_*.nii')
            popd(pwd0);
        end
        function T = build_table_variances(this, subgroups, opts)
            arguments
                this mladni.NMF
                subgroups {mustBeText}
                opts.workdir {mustBeFolder} = fullfile(getenv("ADNI_HOME"), "NMF_FDG")
                opts.numbases double = this.N_PATTERNS
            end
            subgroups = convertCharsToStrings(subgroups);
            niidir = fullfile(opts.workdir, "baseline_cn", "NumBases"+opts.numbases, "OPNMF", "niiImg");

            % NMF patterns
            for ib = 1:opts.numbases
                basis_ics(ib) = mlfourd.ImagingContext2(fullfile(niidir, "Basis_"+ib+".nii")); 
            end

            % table of variances for patterns (rows) and subgroups (cols)
            T = table();
            for isg = 1:length(subgroups)
                % T columns
                var_fns = glob(fullfile(opts.workdir, subgroups(isg), "all_trc-FDG*pet_on_T1w_Warped_dlicv_variance.nii.gz"));
                assert(~isempty(var_fns))
                var_ic = mlfourd.ImagingContext2(var_fns{1});

                var_vec = nan(opts.numbases, 1);
                for ib = 1:opts.numbases
                    % T rows
                    var_vec(ib) = var_ic.volumeWeightedAveraged(basis_ics(ib)./dipsum(basis_ics(ib))); % weighted average
                end

                T = addvars(T, ascol(var_vec), NewVariableNames=subgroups(isg));
            end

            writetable(T, fullfile(opts.workdir, stackstr()+".csv"));
            save(fullfile(opts.workdir, stackstr()+".mat"), "T");
            heatmap(table2array(T));
            ylabel("NMF Patterns");
            xlabel("Diagnostic Groups");
            saveFigure2(gcf, stackstr());
        end
        function call(this)
            %% Writes imaging data to X.mat, 
            %  then calls this.calcRecError(), which plots and saves *RecError.{fig,png}.

            param.isList = 1 ;
            param.downSample = this.downSample;
            param.smooth = 0;
            param.mask = fullfile(getenv('ADNI_HOME'), 'VolBin/mask.nii.gz');
            assert(isfile(param.mask))
            param.permute = 0;
            param.numBase = this.selectedNumBases;
            param.componentDir = this.componentDir;
            
            if isfile(this.Xmat) && this.cache_files
                load(this.Xmat);
            else
                data = this.loadData(this.inFiles, param, []);
                meanX = mean(data.X,2);
                X = data.X(meanX>0,:); 
                clear data;
                save(this.Xmat, 'X', 'meanX');
            end
            resultsDir = fullfile(this.outputDir, 'results');
            ensuredir(resultsDir);

            % BUG FIX:  selectVoxels must be logical, and
            %           numel(selectVoxels) == numel(mask_for_imaging), and
            %           size(selectVoxels) == [N_voxels 1], and
            %           size(X,1) == sum(selectVoxels)
            assert(isfile(param.mask))
            mask_ic = mlfourd.ImagingContext2(param.mask);
            selectVoxels = logical(mask_ic);
            selectVoxels = selectVoxels(:);
            assert(size(X,1) == sum(selectVoxels), stackstr())
            this.calcRecError(X, this.outputDir, resultsDir, selectVoxels);
        end
        function call2(this, tags)
            %% Writes imaging data to X2.mat, 
            %  then calls this.calculateSelectedComponentWeightedAverageNIFTI(), 
            %  which writes component_weighted_average_study-study_design.csv.

            arguments
                this mladni.NMF
                tags {mustBeTextScalar} = this.study_design
            end

            param.isList = 1 ;
            param.downSample = this.downSample;
            param.smooth = 0;
            param.mask = fullfile(getenv('ADNI_HOME'), 'VolBin/mask.nii.gz');
            assert(isfile(param.mask))
            param.permute = 0;
            param.numBase = this.selectedNumBases;
            param.componentDir = this.componentDir;

            if isfile(this.X2mat) && this.cache_files
                load(this.X2mat);
            else
                data = this.loadData(this.inFiles2, param, []);
                meanX = mean(data.X,2);
                X = data.X(meanX>0,:); 
                clear data;
                save(this.X2mat, 'X', 'meanX');
            end

            weightedAverFilename = fullfile(this.componentDir, sprintf('component_weighted_average_%s.csv', tags));
            this.calculateSelectedComponentWeightedAverageNIFTI( ...
                this.inFiles2, this.targetDatasetDir, this.selectedNumBases, weightedAverFilename);
        end
        function run_extractBasesMT(this, numBases, outputDir)
            assert(isscalar(numBases));
            assert(isfolder(outputDir));
            
            cmd = sprintf('%s/run_extractBasesMT.sh %s OPNMF %s 1 %i outputDir %s saveInterm 1 negPos 0 initMeth 4', ...
                this.volbin, this.mcrroot, this.inFiles, numBases, outputDir); % downSample 2
            system(cmd);
        end
        function run_nmf_dataset(this)
            %% quasi-constant

            this.mcrroot = '/data/nil-bluearc/raichle/jjlee/Local/MCR_zfs/R2018b/v95';
            if ~isempty(ipr.numBases)
                this.numBases = ipr.numBases;
            else
                for n = 2:2:mladni.NMF.MAX_NUM_BASES
                    if ~isfile(fullfile(this.outputDir, sprintf('NumBases%i', n), 'OPNMF', 'niiImg', ...
                            sprintf('Basis_%i.nii', n)))
                        try
                            cmd = sprintf('%s/run_nmf_dataset.sh -d %s -r %i', ...
                                this.volbin, this.nmfDataset, n);
                            fprintf(strcat('mladni.NMF: ', cmd))
                            mysystem(cmd);
                        catch ME
                            handerror(ME)
                        end                        
                    end
                    this.numBases = [this.numBases n];
                end
            end
        end
        function submit_nmf_dataset(this)
            cmd = sprintf('%s/submit_nmf_dataset.sh -d %s -m %i', this.volbin, this.nmfDataset, this.memInGB);
            [~,r] = mlbash(cmd);
            disp(r)
        end
        
        function this = NMF(varargin)
            %% NMF works best with pwd ~ getenv('ADNI_HOME')/NMF_FDG
            %  Args:
            %      downSample (double): >= 1
            %      memInGB (scalar):  slurm mem request per batch job.
            %      nmfDataset (text): tag for dataset, used to specify ADNI/NMF_FDG/nmfDataset/nifti_files.csv; 
            %                         default is "test".
            %      volbin (text): home for NMF standalone executable.
            
            ip = inputParser;
            addParameter(ip, "cache_files", false, @islogical);
            addParameter(ip, "downSample", 1, @(x) x >= 1);
            addParameter(ip, "isList", true, @islogical);
            addParameter(ip, "mask", fullfile(getenv("ADNI_HOME"), "VolBin", "mask.nii.gz"), @isfile);
            addParameter(ip, "memInGB", 4, @isscalar);
            addParameter(ip, "nmfDataset", "baseline_cn", @(x) istext(x));
            addParameter(ip, "numBases", 2:2:mladni.NMF.MAX_NUM_BASES, @isnumeric)
            addParameter(ip, "permute", false, @islogical);
            addParameter(ip, "repetitions", 50, @isscalar);
            addParameter(ip, "selectedNumBases", 24, @isscalar);
            addParameter(ip, "smooth", false, @islogical);
            addParameter(ip, "study_design", "longitudinal", @istext)
            addParameter(ip, "volbin", fullfile(getenv("ADNI_HOME"), "VolBin"), @isfolder);
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.cache_files = ipr.cache_files;
            this.downSample = ipr.downSample;
            this.isList = ipr.isList;
            this.mask = ipr.mask;
            this.memInGB = ipr.memInGB;
            this.nmfDataset = ipr.nmfDataset;
            assert(isfolder(this.outputDir));
            this.study_design = ipr.study_design;
            this.permute = ipr.permute;
            this.repetitions = ipr.repetitions;
            this.selectedNumBases = ipr.selectedNumBases;
            this.smooth = ipr.smooth;
            this.volbin = ipr.volbin; 
        end
    end

    %% PRIVATE

    properties (Access = private)
        adniDemographics_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
