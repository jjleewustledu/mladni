 classdef NMF < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% baseline5 ~ 22 bases
    %  baseline4 ~ 28 bases
    %  baseline_pve1_msk ~ 16 bases
    %  baseline_pve1_msk_fdg ~ 14 bases
    %
    %  baseline_cn ~ 16 bases
    %  baseline_preclinical ~ 30 bases
    %  baseline_cdr_0p5_apos_emci ~ 34 bases
    %  baseline_cdr_0p5_apos_lmci ~ 32 bases
    %  baseline_cdr_gt_0p5_apos ~ 38 bases
    %  
    %  Created 23-Jun-2022 13:02:04 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
        
    methods (Static)
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
        function getDebugLog(j,c)
            try
                c.getDebugLog(j)
            catch
                c.getDebugLog(j.Tasks(end))
            end
        end
        function create()
            dx = {'cn', 'preclinical', 'cdr_0p5_aneg', 'cdr_0p5_apos', 'cdr_0p5_anan', 'cdr_gt_0p5_apos'};
            rep = {'repA', 'repB'};
            parfor m = 2:length(dx)
                for n = 1:2
                    set = sprintf('baseline_%s_%s', dx{m}, rep{n});
                    pth = fullfile(getenv('ADNI_HOME'), 'NMF_FDG', set);
                    pwd0 = pushd(pth);
                    this = mladni.NMF(nmfDataset=set);
                    call(this)
                    popd(pwd0);
                end
            end
        end
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
            if strcmpi(basename(ipr.path), 'niiImg') % /path/to/niiImg
                pwd0 = pushd(ipr.path);
                atl = mlfourd.ImagingContext2( ...
                    fullfile(getenv('FSLDIR'), 'data', 'standard', 'MNI152_T1_1mm.nii.gz'));
                for g = glob('Basis_*.nii')'
                    ic = mlfourd.ImagingContext2(g{1});
                    product = strrep(strcat(ic.fqfp, '.png'), 'niiImg', 'Figures');
                    if ~isfile(product) || ...
                            (isfile(product) && ~ipr.noclobber)
                        atl.save_qc(ic);
                        png = strrep(g{1}, '.nii', '.png');
                        movefile(png, '../Figures', 'f');
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
            ip = inputParser;
            addRequired(ip, 'csv', @isfile);
            addParameter(ip, 'filepath', pwd, @isfolder);
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
            ic.filepath = ipr.filepath;
            ic.save();
            
            % create variance
            ic1 = mlfourd.ImagingContext2(char(csv(1)));
            ic1.selectMatlabTool();
            ic1 = (ic1 - ic).^2;
            re1 = regexp(ic.fileprefix, '(?<descrip>\S+)_mean', 'names');
            fp1 = sprintf('all_%s_variance', re1.descrip);
            for i = 2:length(csv)
                try
                    ic1 = ic1 + (mlfourd.ImagingContext2(char(csv(i))) - ic).^2;
                    ic1.fileprefix = fp1;
                catch
                end
            end
            ic1 = ic1/(length(csv) - 1);
            ic1.fileprefix = fp1;
            ic1.filepath = ipr.filepath;
            ic1.save();
        end
        function create_median_imaging(varargin)
            ip = inputParser;
            addRequired(ip, 'csv', @isfile);
            addParameter(ip, 'filepath', pwd, @isfolder);
            parse(ip, varargin{:});
            ipr = ip.Results;            
            
            csv = readlines(ipr.csv); % string vector
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
            ifc.filepath = ipr.filepath;
            ifc.save();            
        end
        function csv_out = filter_missing_images(varargin)
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
                    lines1 = [lines1; lines{li}];
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
        function parcall()
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
            targetDatasetDir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG', 'baseline_cn');
            assert(isfolder(targetDatasetDir))

            b = 2:2:40;
            parfor bi = 1:length(b)
                groups = mladni.NMF.groups1;
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
                parfor gi = 1:length(groups)
                    this = mladni.NMF(nmfDataset=sprintf('baseline_%s', groups{gi}));
                    cd(this.outputDir);
                    this.create_montage(this.outputDir, noclobber=true);
                end
            else
                groups = mladni.NMF.groups0;
                parfor gi = 1:length(groups)
                    this = mladni.NMF(nmfDataset=sprintf('baseline_%s', groups{gi}), ...
                        workDir=fullfile(getenv("ADNI_HOME"), "NMF_FDG", "Previous"));
                    cd(this.outputDir);
                    this.create_montage(this.outputDir, noclobber=true);
                end
            end
        end
        function split_nifti_files(csv)
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
        
        %% Aris' inference
        
        function calculateComponentWeightedAverageNIFTI(dataList,resultsDir,numBases,outPath)
                        
            % dataList: .csv file containing the images for which the coefficients need
            %           to be calculated. Full path for every image file is given in
            %           every line
            % resultsDir: directory where the NMF results are (i.e., the level where
            %             the NumBases folder is placed)
            % numBases: determines the solution for which one wants to calculate
            %           subject coefficients
            % outPath: CSV path to save output
            
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
        function calculateRepeatedComponentWeightedAverageNIFTI(dataLists,resultsDir,numBases,outFile)
                        
            % dataList: .csv file containing the images for which the coefficients need
            %           to be calculated. Full path for every image file is given in
            %           every line
            % resultsDir: directory where the NMF results are (i.e., the level where
            %             the NumBases folder is placed)
            % numBases: determines the solution for which one wants to calculate
            %           subject coefficients
            % outPath: CSV path to save output
            
            import mlniftitools.*;
            resultsDir = convertStringsToChars(resultsDir);
            outFile = convertStringsToChars(outFile);
        
            disp(numBases)
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
                    
                    %print statement added by tom
                    disp(sprintf('(%d/%d)', i, count)); %#ok<DSPS>
                    
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
        function RecError = calcRecError(X, resultDir, outputDir, selectVoxels)
            
            % Function that calculates the reconstruction error given data matrix X and the non-negative matrix 
            % factorization results saved in the resultDir directory. The function additionally saves figures that plot
            % the reconstruction error, the gradient of the reconstruction error, and the percentage of improvement as a 
            % function of the number of components
            
            % load results and calculate reconstruction error
            % We assume that the results directory is organized as follows: there is a folder for every solution (i.e.,
            % NumBases${K}, wehere K is the number of components for the solution), and inside the folder there is .mat 
            % file (ResultsExtractBases.mat) that contains the matrices W and H that were estimated by the non negative 
            % matrix factorization

            resultDir = convertStringsToChars(resultDir);
            outputDir = convertStringsToChars(outputDir);
            
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
                load([resultDir '/NumBases' num2str(sortedBasisNum(b)) '/OPNMF/ResultsExtractBases.mat']) %#ok<LOAD>
                Est = B*C ;
                Est = Est(selectVoxels,:);
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
        function evaluateRepeatedReproducibility(subgroup, N, K, saveGcf)
            arguments
                subgroup {mustBeTextScalar} = 'cn'
                N double {mustBeInteger} = 20 % # repetitions
                K double {mustBeInteger} = 20 % # bases checked by VolBin/*
                saveGcf logical = true
            end

            home = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
            subgroups = { ...
                'cn', 'preclinical', ...
                'cdr_0p5_aneg_emci', 'cdr_0p5_aneg_mci', 'cdr_0p5_aneg_lmci', ...
                'cdr_0p5_apos_emci', 'cdr_0p5_apos_mci', 'cdr_0p5_apos_lmci', ...
                'cdr_gt_0p5_aneg', 'cdr_gt_0p5_apos'};
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
                N double {mustBeInteger} = 20 % # repetitions
                K double {mustBeInteger} = 20 % # bases checked by VolBin/*
                saveGcf logical = true
                opts.fracx double = 0.5
                opts.Npx double = 3400
                opts.fracy double = 0.5
                opts.Npy double = 1440
            end

            import mladni.Jones2022.boxchart

            home = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
            subgroups = { ...
                'cn', 'preclinical', ...
                'cdr_0p5_aneg_emci', 'cdr_0p5_aneg_mci', 'cdr_0p5_aneg_lmci', ...
                'cdr_0p5_apos_emci', 'cdr_0p5_apos_mci', 'cdr_0p5_apos_lmci', ...
                'cdr_gt_0p5_aneg', 'cdr_gt_0p5_apos'};
            assert(any(strcmp(subgroup, subgroups)));
            
            outputDir_ = fullfile(home, sprintf('baseline_%s', subgroup));
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
                            catch
                            end
                        end
                    catch ME
                        handwarning(ME)
                    end
                end
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
                N_ARI = numel(ARI);                
                basisLabel = repmat(sortedBasisNum', [1 N]);                

                figure
                boxchart(reshape(basisLabel, [N_ARI 1]), sortedBasisNames, reshape(ARI, [N_ARI 1]), 4);
                xlabel('Number of components','fontsize',20)
                ylabel({'Split-sample reproducibility' ;'(Adjusted Rand Index)'},'fontsize',20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf,[outputDir_ '/AdjustedRandIndexReproducibility_repeat2.fig'])
                saveas(gcf,[outputDir_ '/AdjustedRandIndexReproducibility_repeat2.png'])

                figure
                boxchart(reshape(basisLabel, [N_ARI 1]), sortedBasisNames, reshape(meanInner, [N_ARI 1]), 4);
                xlabel('Number of components','fontsize',20)
                ylabel({'Split-sample reproducibility';'(mean inner product)'},'fontsize',20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf,[outputDir_ '/MeanInnerProductReproducibility_repeat2.fig'])
                saveas(gcf,[outputDir_ '/MeanInnerProductReproducibility_repeat2.png'])

                figure
                boxchart(reshape(basisLabel, [N_ARI 1]), sortedBasisNames, reshape(medianInner, [N_ARI 1]), 4);
                xlabel('Number of components','fontsize',20)
                ylabel({'Split-sample reproducibility';'(median inner product)'},'fontsize',20)
                set(gca,'fontsize',20)
                set(gcf, Position=[1 1 opts.fracx*opts.Npx opts.fracy*opts.Npy])
                saveas(gcf,[outputDir_ '/MedianInnerProductReproducibility_repeat2.fig'])
                saveas(gcf,[outputDir_ '/MedianInnerProductReproducibility_repeat2.png'])

            end
        end
        function [meanInner,medianInner,ARI,overlap,sortedBasisNum] = evaluateReproducibility( ...
                pathDir1, pathDir2, outputDir, saveGcf)
            %% https://github.com/sotiraslab/aris_nmf_analyses/blob/main/evaluateReproducibility.m
            %
            % function that evaluates the reproducibility between the results obtained
            % for two different splits
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
                data = loadDataFromList(dataInput,param,subsetIdx) ; % data is a struct !
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
    
    properties
        inFiles
        mcrroot
        memInGB
        nmfDataset
        numBases
        selectedNumBases
        volbin
        Xmat
        X2mat
    end

    properties (Constant)
        groups0 = { ...
            'cdr_0p5_apos_emci', 'cdr_0p5_apos_lmci'}
        groups = { ...
            'cn', 'preclinical', ...
            'cdr_gt_0p5_apos', ...
            'cdr_0p5_apos_emci', 'cdr_0p5_apos_lmci', ...
            'cdr_ge_0p5_aneg', ... 
            'cdr_0p5_aneg_emci', 'cdr_0p5_aneg_lmci', ...
            'cdr_0p5_aneg', 'cdr_0p5_apos', 'cdr_0p5_anan'}
        groups1 = { ...
            'cn', 'preclinical', ...
            'cdr_0p5_apos_emci', 'cdr_0p5_apos_lmci', ...
            'cdr_gt_0p5_apos', ...
            'cdr_ge_0p5_aneg'}
        bases1 = [16 30 34 32 38]
        repetitions = 20
    end

    properties (Dependent)
        outputDir
    end

    methods % GET
        function g = get.outputDir(this)
            g = fullfile(getenv('ADNI_HOME'), 'NMF_FDG', this.nmfDataset, '');
        end
    end
    
    methods
        function call(this)
            componentDir = fullfile(this.outputDir, sprintf('NumBases%i', this.selectedNumBases), 'components');
            ensuredir(componentDir);

            param.isList = 1 ;
            param.downSample = 1 ;
            param.smooth = 0;
            param.mask = [];
            param.permute = 0;
            param.numBase = this.selectedNumBases;
            param.componentDir = componentDir;

            if isfile(this.Xmat)
                load(this.Xmat);
            else
                data = this.loadData(this.inFiles,param,[]);
                meanX = mean(data.X,2);
                X = data.X(meanX>0,:); 
                clear data;
                save(this.Xmat, 'X', 'meanX');
            end
            resultsDir = fullfile(this.outputDir, 'results');
            ensuredir(resultsDir);
            this.calcRecError(X, this.outputDir, resultsDir, meanX>0);

            this.calculateComponentWeightedAverageNIFTI( ...
                this.inFiles, this.outputDir, this.selectedNumBases, fullfile(componentDir, 'component_weighted_average.csv'));

            pwd0 = pushd(componentDir);
            fdg = mladni.FDGDemographics();
            fdg.table_covariates('component_weighted_average.csv');
            popd(pwd0);
        end
        function call2(this, targetDatasetDir, tags)
            arguments
                this mladni.NMF
                targetDatasetDir {mustBeFolder} = this.outputDir
                tags {mustBeTextScalar} = ''
            end

            param.isList = 1 ;
            param.downSample = 1 ;
            param.smooth = 0;
            param.mask = fullfile(getenv('ADNI_HOME'), 'VolBin/mask.nii.gz');
            assert(isfile(param.mask))
            param.permute = 0;
            param.numBase = this.selectedNumBases;

            % cache file
            if isfile(this.X2mat)
                load(this.X2mat);
            else
                data = this.loadData(this.inFiles,param,[]);
                meanX = mean(data.X,2);
                X = data.X(meanX>0,:); 
                clear data;
                save(this.X2mat, 'X', 'meanX');
            end

            componentDir = fullfile(this.outputDir, sprintf('NumBases%i', param.numBase), 'components');
            ensuredir(componentDir);
            csvFilename = sprintf('component_weighted_average_%s.csv', tags);
            this.calculateRepeatedComponentWeightedAverageNIFTI( ...
                this.inFiles, targetDatasetDir, this.selectedNumBases, fullfile(componentDir, csvFilename));

            pwd0 = pushd(componentDir);
            fdgd = mladni.FDGDemographics();
            fdgd.table_covariates(csvFilename, tags=tags);
            popd(pwd0);
        end
        function inFiles2 = replace_inFiles(~, inFiles)
            
            tbl = readtable(inFiles, 'ReadVariableNames', false, 'Delimiter', ' ');
            Var1 = tbl.Var1;
            Var1 = strrep(Var1, '/home/aris_data/ADNI_FDG', getenv('ADNI_HOME'));
            Var1 = strrep(Var1, '/scratch/jjlee/Singularity/ADNI', getenv('ADNI_HOME'));
            
            [pth,fp,x] = myfileparts(inFiles);
            inFiles2 = fullfile(pth, strcat(fp, '2', x));
            writetable(table(Var1), inFiles2, 'WriteVariableNames', false);
        end
        function run_extractBasesMT(this, numBases, outputDir)
            assert(isscalar(numBases));
            assert(isfolder(outputDir));
            
            cmd = sprintf('%s/run_extractBasesMT.sh %s OPNMF %s 1 %i outputDir %s saveInterm 1 negPos 0 initMeth 4', ...
                this.volbin, this.mcrroot, this.inFiles, numBases, outputDir); % downSample 2
            system(cmd);
        end
        function submit_nmf_dataset(this)
            cmd = sprintf('%s/submit_nmf_dataset.sh -d %s -m %i', this.volbin, this.nmfDataset, this.memInGB);
            [~,r] = mlbash(cmd);
            disp(r)
        end
        
        function this = NMF(varargin)
            %% NMF works in getenv('SINGULARITY_HOME')
            %  Args:
            %      memInGB (scalar):  slurm mem request per batch job.
            %      nmfDataset (text): tag for dataset, used to specify workDir/nmfDataset/nifti_files.csv; 
            %                         default is "test".
            %      workDir (text): containing subfolder for dataset, then subfolders for variable NMF components.
            %      volbin (text): home for NMF standalone executable.
            
            if contains(hostname, 'pascal')
                setenv('SINGULARITY_HOME', '/mnt/CHPC_scratch/Singularity');
            end
            if contains(hostname, 'vglab')
                setenv('SINGULARITY_HOME', '/home/usr/jjlee/mnt/CHPC_scratch/Singularity');
            end

            ip = inputParser;
            addParameter(ip, "memInGB", 4, @isscalar);
            addParameter(ip, "nmfDataset", "baseline_cn", @(x) istext(x));
            addParameter(ip, "numBases", 2:2:40, @isnumeric)
            addParameter(ip, "selectedNumBases", 18, @isscalar);
            addParameter(ip, "volbin", ...
                fullfile(getenv("ADNI_HOME"), "VolBin"), @isfolder);
            addParameter(ip, "workDir", ...
                fullfile(getenv("ADNI_HOME"), "NMF_FDG"), @isfolder);
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.memInGB = ipr.memInGB;
            this.nmfDataset = ipr.nmfDataset;
            assert(isfolder(this.outputDir));
            this.inFiles = fullfile(this.outputDir, 'nifti_files.csv');
            this.inFiles = this.replace_inFiles(this.inFiles);
            this.inFiles = this.filter_missing_images(this.inFiles);
            assert(isfile(this.inFiles));
            this.selectedNumBases = ipr.selectedNumBases;
            this.volbin = ipr.volbin;
            numBasesFolder = sprintf('NumBases%i', ipr.selectedNumBases);
            this.Xmat = fullfile(this.outputDir, numBasesFolder, 'X.mat');
            this.X2mat = fullfile(this.outputDir, numBasesFolder, 'X2.mat');
            
            % quasi-constant
            this.mcrroot = '/data/nil-bluearc/raichle/jjlee/Local/MCR_zfs/R2018b/v95';
            if ~isempty(ipr.numBases)
                this.numBases = ipr.numBases;
            else
                for n = 2:2:40
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
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
