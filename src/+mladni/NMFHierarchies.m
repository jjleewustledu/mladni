classdef NMFHierarchies < handle
    %% line1
    %
    %  See also this.table_patt_weighted_fdg for ranking by component-weighted averages of baseline_cn.
    %  
    %  Created 06-Feb-2024 20:51:19 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Constant)
        EPS = 1e-6
        N_PATTERNS = mladni.NMF.N_PATTERNS
        selected_spans = [2, 6, 8, 10, 12, 14, 24]
    end

    properties    
        data_home
        home
        matfile0 = 'NMFCovariates_table_covariates_1stscan_longitudinal.mat'
        Nforeground = 228483  % mask has 228483 foreground voxels in 2mm, 67019 voxels in 3mm
    end

    properties (Dependent)
        mapping_span_by_suvr
        mask
        standard_dir
        workdir
    end

    methods %% GET
        function g = get.mapping_span_by_suvr(this)
            if ~isempty(this.mapping_span_by_suvr_)
                g = this.mapping_span_by_suvr_;
                return
            end

            T = this.table_patt_weighted_fdg();
            T = sortrows(T, "indices_bases");
            this.mapping_span_by_suvr_ = T.rank;
            g = this.mapping_span_by_suvr_;
        end
        function g = get.mask(this)
            mask_fqfn = fullfile(getenv("ADNI_HOME"), "VolBin", "mask.nii.gz");
            assert(isfile(mask_fqfn))
            g = mlfourd.ImagingContext2(mask_fqfn);
            assert(sum(g.imagingFormat.img, "all") == this.Nforeground)
        end
        function g = get.standard_dir(~)
            g = fullfile(getenv('FSLDIR'), 'data', 'standard');
        end
        function g = get.workdir(this)
           g = fullfile(this.data_home, 'NMF_FDG');
        end
    end

    methods
        function this = NMFHierarchies(opts)
            %% See also this.table_patt_weighted_fdg for ranking by component-weighted averages of baseline_cn.

            arguments
                opts.data_home {mustBeFolder} = getenv("ADNI_HOME")
                opts.home {mustBeTextScalar} = ""  % e.g. baseline_cn
            end

            this.data_home = opts.data_home;
            this.home = opts.home;
            if isemptytext(opts.home)
                this.home = fullfile(this.workdir, "baseline_cn");
            end
        end

        function build_argmax_maps(this, varargin)
            %% builds argmax maps for each of the spanning spaces for a cohort, e.g., baseline_cn.

            pwd0 = pushd(this.home);
            for s = 2:2:max(this.selected_spans)
                pwd1 = pushd(fullfile("NumBases"+s, "OPNMF", "niiImg")); 
                this.build_argmax_map(s, varargin{:}); 
                popd(pwd1); 
            end
            popd(pwd0);
        end

        function build_argmax_map(this, span, opts)
            %% Identifies voxel intensities with argmax, arg labeling the basis;
            %  also sorts args with arg==1 having greatest total SUVR.
            %  See also this.table_patt_weighted_fdg for ranking by component-weighted averages of baseline_cn.

            arguments
                this mladni.NMFHierarchies
                span {mustBeNumeric}
                opts.rank_by_suvr logical = true
            end
            
            ic = mlfourd.ImagingContext2('Basis_1.nii');
            ic = ic.zeros;
            ifc = ic.nifti;
            for b = 1:span
                ifc1 = mlfourd.ImagingFormatContext2(sprintf('Basis_%i.nii', b));
                if opts.rank_by_suvr
                    b1 = this.ranking_by_suvr(span, b);
                else
                    b1 = b;
                end
                ifc.img(:,:,:,b1) = ifc1.img;
            end
            ifc.fileprefix = 'Basis_all'; % composite nifti of bases
            ifc.save

            mat = reshape(ifc.img, [91*109*91, span]);  % N_voxels x K
            mat = this.normalize_cols_rows(mat);  % cols <- cols/norm(col); rows <- rows/sum(rows)
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

        function [Ts, Us] = build_component_weighted_averages_and_covariates(this)
            Ts = cell(1, length(this.selected_spans));
            Us = cell(1, length(this.selected_spans));
            for idx = 1:length(this.selected_spans)-1
                fprintf("%s: this.selected_spans(%i) -> %i\n", stackstr(), idx, this.selected_spans(idx))
                nmf = mladni.NMF(selectedNumBases=this.selected_spans(idx));
                Ts{idx} = nmf.call2();
                Us{idx} = this.table_covariates_1stscan(selectedNumBases=this.selected_spans(idx));
            end                     
        end

        function ic = build_downsampled_to_3mm(this, ic)
            arguments
                this mladni.NMFHierarchies
                ic mlfourd.ImagingContext2
            end

            refdir = getenv("REFDIR");
            ref = mlfourd.ImagingContext2(fullfile(refdir, "MNI152_T1_3mm.nii.gz"));
            out = mlfourd.ImagingContext2(fullfile(ic.fqfp + "_3mm.nii.gz"));
            omat = fullfile(fullfile(refdir, "MNI152_T1_2mm_on_3mm.mat"));
            flirt = mlfsl.Flirt( ...
                in=ic, ...
                ref=ref, ...
                out=out, ...
                omat=omat, ...
                interp='nearestneighbour', ...
                noclobber=false);
            %flirt.flirt();
            flirt.applyXfm();
            assert(isfile(out.fqfn));
            ic = out;
        end

        function T = build_table_for_ggalluvial2(this)
            %% maintain voxel identities across pattern-models

            pwd0 = pushd(this.home);            

            Nspans = length(this.selected_spans);
            A = NaN(this.Nforeground, Nspans);
            for sidx = length(this.selected_spans):-1:1
                mg = mglob(sprintf("NumBases%i/OPNMF/niiImg/Basis_argmax_brain_mask.nii", this.selected_spans(sidx)));
                ic = mlfourd.ImagingContext2(mg(1));
                img = ic.imagingFormat.img;
                img = ascol(img(img > 0));
                assert(length(img) == this.Nforeground)
                
                if sidx == length(this.selected_spans)
                    img_span_ref = img;
                else
                    img = this.reorder_patterns(img, img_span_ref);
                end

                ifc = copy(this.mask.imagingFormat);
                ifc.img(ifc.img > 0) = img;
                ifc.fqfileprefix = ic.fqfileprefix + "_reordered_patterns";
                ifc.save();

                A(:, sidx) = img;  % A <- assignment of pattern-# in pattern-model to voxel
            end

            tag_ = A(:, end);  % tag_ <- span-24
            % tag = cell(size(tag_));
            % for ti = 1:length(tag_)
            %     tag{ti} = sprintf('P%02i', tag_(ti));
            % end
            tag = string(tag_);

            freq = ones(this.Nforeground, 1);

            T = table(A, tag, freq, VariableNames=["A", "Tag", "freq"]);
            T = splitvars(T, "A", NewVariableNames="N_Patterns_"+asrow(string(this.selected_spans)));
            writetable(T, stackstr()+".csv")

            popd(pwd0)
        end

        function img1 = reorder_patterns(this, img, ref)
            %% Use cosine similarity to reorder pattern indices in img to have best similarity to img_span_ref, e.g. ,
            %  NumBasis24.

            arguments
                this mladni.NMFHierarchies
                img {mustBeNumeric}
                ref {mustBeNumeric}
            end

            if all(img == ref)
                return
            end
            img = asrow(img);
            ref = asrow(ref);
            span_img = max(img, [], "all");
            span_ref = max(ref, [], "all");

            img1 = img;
            sim = nan(1, span_ref);
            for q = 1:span_img
                for p = 1:span_ref
                    pref = double(ref == p);
                    qimg = double(img == q);
                    sim(p) = this.cosine_similarity([pref; qimg]);
                end
                [~,new_idx] = max(sim);
                img1(img == q) = new_idx;
            end            
        end

        function sim = cosine_similarity(~, A)
            sim = 1 - pdist(A, 'cosine');
        end

        function mu = mu_suvr_for_basis(this, b, opts)
            arguments
                this mladni.NMFHierarchies
                b double
                opts.span double = 24
            end

            T = this.table_patt_weighted_fdg(span=opts.span);
            mu = T.mu(T.indices_bases == b);
        end

        function g = ranking_by_suvr(this, span, basis)
            arguments
                this mladni.NMFHierarchies
                span double
                basis double
            end

            if ~isempty(this.ranking_by_suvr_)
                T = this.ranking_by_suvr_.("table_span"+span);
                g = T{T.indices_bases == basis, 'rank'};
                return
            end

            this.ranking_by_suvr_ = struct();
            for ssi = 2:2:max(this.selected_spans)
                span_ = ssi;
                T_ = this.table_patt_weighted_fdg(span=span_);
                this.ranking_by_suvr_.("table_span"+span_) = T_;
            end

            T = this.ranking_by_suvr_.("table_span"+span);
            g = T{T.indices_bases == basis, 'rank'};
        end

        function T = table_covariates_1stscan(this, opts)
            arguments
                this mladni.NMFHierarchies
                opts.selectedNumBases {mustBeScalarOrEmpty} = this.N_PATTERNS
            end
            nmfc = mladni.NMFCovariates(selectedNumBases=opts.selectedNumBases);            

            % also save NMFCovariates.covariates_1stscan_file
            T = nmfc.table_covariates_1stscan();
        end

        function T = table_patt_weighted_fdg(this, opts)

            % nmfh.table_patt_weighted_fdg
            % ans =
            %   24×3 table
            %
            %  rank  indices_bases    mu        sigma
            %  ____  _____________  _______    ________
            %
            %  1        22           1.3149     0.12809
            %  2        10           1.2628     0.13907
            %  3        11           1.2406     0.15168
            %  4        15           1.2227     0.13558
            %  5         5            1.206     0.12444
            %  6        17           1.2051     0.12955
            %  7        13           1.1859     0.14648
            %  8         9           1.1639     0.13782
            %  9         3           1.1386     0.14355
            %  10       24           1.1309     0.11426
            %  11       16           1.1286     0.10461
            %  12        2           1.1146     0.12517
            %  13        4           1.0951     0.10383
            %  14        8           1.0474    0.075412
            %  15       14           1.0423     0.10174
            %  16       20           1.0375    0.036501
            %  17        7           1.0082    0.099031
            %  18       12           1.0022     0.12901
            %  19       18          0.99447     0.15053
            %  20       23          0.96043    0.067123
            %  21        1          0.95748    0.097545
            %  22       19          0.91651    0.067689
            %  23       21          0.77652     0.10514
            %  24        6          0.72945    0.060317
            
            arguments
                this mladni.NMFHierarchies
                opts.span = this.N_PATTERNS
            end
            assert(any(ismember(opts.span, 1:max(this.selected_spans))))

            c = 1;
            ld = load(fullfile( ...
                this.workdir, ...
                'baseline_cn', ...
                sprintf('NumBases%i', opts.span), ...
                'components', ...
                this.matfile0));

            indices_bases = 1:opts.span; 
            mu = nan(1, opts.span);
            sigma = nan(1, opts.span);
            for idx = 1:opts.span
                comp = ld.t.(sprintf("Components_%i", idx));
                mu(idx) = mean(comp, 1);
                sigma(idx) = std(comp, 1);
            end

            T = table(indices_bases', mu', sigma', VariableNames={'indices_bases', 'mu', 'sigma'}); 
            T = sortrows(T, 2, "descend");
            T = addvars(T, ascol(1:length(indices_bases)), NewVariableNames={'rank'});
            this.sorted_bases_ = T.indices_bases;
        end
    end

    methods (Static)
        function M = normalize_cols_rows(M)
            %% prescribed per discussion with Aris 4/15/2024.
            %  M ~ 1e6 x 24, so normalize only nonzero voxels.

            nonzero = sum(M, 2) > 0;

            % cols <- cols/norm(col)
            for n = 1:size(M, 2)
                V = M(nonzero, n);
                V = V / norm(V);
                M(nonzero, n) = V;
            end
            
            % rows <- rows/sum(rows); 
            M1 = M(nonzero, :);
            for m = 1:size(M1, 1)
                S = M1(m, :);
                S = S / sum(S);
                M1(m, :) = S;
            end
            M(nonzero, :) = M1;
        end
    end

    %% PRIVATE

    properties (Access = private)
        mapping_span_by_suvr_
        ranking_by_suvr_
        sorted_bases_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
