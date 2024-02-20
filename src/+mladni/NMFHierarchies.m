classdef NMFHierarchies < handle
    %% line1
    %
    %  See also NMFRadar.table_patt_weighted_fdg for ranking by component-weighted averages of baseline_cn.
    %  
    %  Created 06-Feb-2024 20:51:19 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Constant)
        EPS = 1e-6
        N_PATTERNS = mladni.NMF.N_PATTERNS
    end

    properties        
        home
        Nforeground = 228483  % mask has 228483 foreground voxels
        selected_spans = [2, 8, 10, 12, 14 24]
    end

    properties (Dependent)
        ranking_by_suvr  % See also NMFRadar.table_patt_weighted_fdg for ranking by component-weighted averages of baseline_cn.
        standard_dir
        mu_suvr_24
    end

    methods %% GET
        function g = get.ranking_by_suvr(this)
            if ~isempty(this.ranking_by_suvr__)
                g = this.ranking_by_suvr__;
                return
            end

            m = containers.Map(KeyType="double", ValueType="double"); % OPNMF index -> suvr ranking by component-weighted averaging of baseline_cn
            m(22) = 1; 
            m(10) = 2;
            m(11) = 3;
            m(15) = 4;
            m(5) = 5;
            m(17) = 6;
            m(13) = 7;
            m(9) = 8;
            m(3) = 9;
            m(24) = 10;
            m(16) = 11;
            m(2) = 12;
            m(4) = 13;
            m(8) = 14;
            m(14) = 15;
            m(20) = 16;
            m(7) = 17;
            m(12) = 18;
            m(18) = 19;
            m(23) = 20;
            m(1) = 21;
            m(19) = 22;
            m(21) = 23;
            m(6) = 24;
            this.ranking_by_suvr__ = m;
            g = m;
        end
        function g = get.standard_dir(~)
            g = fullfile(getenv('FSLDIR'), 'data', 'standard');
        end
        function g = get.mu_suvr_24(this)
            if ~isempty(this.mu_suvr_24__)
                g = this.mu_suvr_24__;
                return
            end

            nmfr = mladni.NMFRadar;
            this.mu_suvr_24__ = nmfr.table_patt_weighted_fdg.mu;
            g = this.mu_suvr_24__;
        end
    end

    methods
        function build_argmax_maps(this)
            %% builds argmax maps for each of the spanning spaces for a cohort, e.g., baseline_cn.

            pwd0 = pushd(this.home);
            for s = 2:2:40
                pwd1 = pushd(fullfile("NumBases"+s, "OPNMF", "niiImg")); 
                this.build_argmax_map(s); 
                popd(pwd1); 
            end
            popd(pwd0);
        end
        function build_argmax_map(this, span)
            %% Identifies voxel intensities with argmax, arg labeling the basis;
            %  also sorts args with arg==1 having greatest total SUVR.
            %  See also NMFRadar.table_patt_weighted_fdg for ranking by component-weighted averages of baseline_cn.
            
            ic = mlfourd.ImagingContext2('Basis_1.nii');
            ic = ic.zeros;
            ifc = ic.nifti;
            if span == this.N_PATTERNS
                for b = 1:span
                    ifc1 = mlfourd.ImagingFormatContext2(sprintf('Basis_%i.nii', b));
                    b1 = this.ranking_by_suvr(b);
                    ifc.img(:,:,:,b1) = ifc1.img;
                end
            else
                for b = 1:span
                    ifc1 = mlfourd.ImagingFormatContext2(sprintf('Basis_%i.nii', b));
                    ifc.img(:,:,:,b) = ifc1.img;
                end
            end
            ifc.fileprefix = 'Basis_all'; % composite nifti of bases
            ifc.save

            mat = reshape(ifc.img, [91*109*91, span]);
            [~,argmax] = max(mat'); %#ok<UDIM> 
            argmax = reshape(argmax, [91, 109, 91]);
            ifc_argmax = copy(ifc);
            ifc_argmax.img = argmax;
            ifc_argmax.fileprefix = 'Basis_argmax';
            ifc_argmax.save();
            
            brain_mask = mlfourd.ImagingContext2(fullfile(this.standard_dir, 'MNI152_T1_2mm_brain_mask.nii.gz'));  % tight mask
            ic_argmax = mlfourd.ImagingContext2(ifc_argmax);
            ic_argmax = ic_argmax .* brain_mask.binarized();
            ic_argmax.fileprefix = strcat(ifc_argmax.fileprefix, '_brain_mask');
            ic_argmax.save();            
        end
        function Ts = build_component_weighted_averages(this)
            Ts = cell(1, length(this.selected_spans));
            for idx = 1:length(this.selected_spans)-1
                fprintf("%s: this.selected_spans(%i) -> %i\n", stackstr(), idx, this.selected_spans(idx))
                nmf = mladni.NMF(selectedNumBases=this.selected_spans(idx));
                Ts{idx} = nmf.call2();
            end                     
        end
        function T = build_table_for_ggalluvial(this)
            pwd0 = pushd(this.home);            

            Nspans = length(this.selected_spans);
            A = NaN(this.Nforeground, Nspans);
            for sidx = 1:length(this.selected_spans)
                mg = mglob(sprintf("NumBases%i/OPNMF/niiImg/Basis_argmax_brain_mask.nii", this.selected_spans(sidx)));
                ifc = mlfourd.ImagingFormatContext2(mg(1));
                img = ifc.img;
                img = ascol(img(img > 0));
                assert(length(img) == this.Nforeground)
                A(:, sidx) = img;
            end

            tag = string(A(:, end));
            tag(tag == "16") = "P16";
            tag(tag == "3") = "P3";
            tag(tag == "2") = "P2";
            tag(tag == "7") = "P7";
            tag(~contains(tag, "P")) = "P";

            suvr = A(:, end);
            for b = 1:this.N_PATTERNS
                suvr(suvr == b) = this.mu_suvr_24(b);
            end

            freq = ones(this.Nforeground, 1);

            T = table(A, tag, suvr, freq, VariableNames=["A", "Tag", "suvr", "freq"]);
            T = splitvars(T, "A", NewVariableNames="Span-"+asrow(string(this.selected_spans)));
            writetable(T, stackstr()+".csv")

            popd(pwd0)
        end

        function this = NMFHierarchies(opts)
            %% See also NMFRadar.table_patt_weighted_fdg for ranking by component-weighted averages of baseline_cn.

            arguments
                opts.home {mustBeFolder} = pwd  % e.g. baseline_cn
            end

            this.home = opts.home;
        end
    end

    methods (Static)
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
                
                % % % find correspondences
                % % [Matching,~] = Hungarian(dist);
                % % [~,idx_hug1]=max(Matching,[],2);
                % % 
                % % % overlap - hungarian
                % % overlap{exp_} = zeros(length(wlen1),1) ;
                % % for b=1:length(wlen1)
                % %     overlap{exp_}(b) = inner_product(b,idx_hug1(b));
                % % end
                
                % % % overlap with best
                % % overlap_best{exp_} = max(inner_product,[],2) ;
                
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
    end

    %% PRIVATE

    properties (Access = private)
        ranking_by_suvr__
        mu_suvr_24__
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
