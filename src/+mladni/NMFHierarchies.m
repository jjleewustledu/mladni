classdef NMFHierarchies < handle
    %% line1
    %
    %  See also this.table_patt_weighted_fdg for ranking by component-weighted averages of baseline_cn.
    %  
    %  Created 06-Feb-2024 20:51:19 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 23.2.0.2485118 (R2023b) Update 6 for MACA64.  Copyright 2024 John J. Lee.
    

    properties (Constant)
        EPS = 1e-6
    end

    properties    
        data_home
        home
        matfile0 = 'NMFCovariates_table_cn_1stscan_longitudinal.mat'
        N_patterns
        Nforeground = 228483  % mask has 228483 foreground voxels in 2mm, 67019 voxels in 3mm
        selected_spans = [2, 6, 8, 10, 12, 14, 24]
    end

    properties (Dependent)
        mask
        standard_dir
        workdir
    end

    methods %% GET
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
                opts.N_patterns double = mladni.NMF.N_PATTERNS
            end

            this.N_patterns = opts.N_patterns;
            this.data_home = opts.data_home;
            this.home = opts.home;
            if isemptytext(opts.home)
                this.home = fullfile(this.workdir, "baseline_cn");
            end
        end

        function build_argmax_maps(this, varargin)
            %% builds argmax maps for each of the spanning spaces for a cohort, e.g., baseline_cn.

            pwd0 = pushd(this.home);
            for s = max(this.selected_spans):-2:2

                nmfc = mladni.NMFCovariates(N_patterns=s);
                t = nmfc.table_cn(true);  % cross-sectional=true
                pwd1 = pushd(fullfile("NumBases"+s, "components"));
                save("NMFCovariates_table_cn_1stscan_longitudinal.mat", "t");
                writetable(t, "NMFCovariates_table_cn_1stscan_longitudinal.csv", WriteVariableNames=true)
                popd(pwd1);

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
                    b1 = this.pattern_num_from(span, b);
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

        function build_maps_for_ggalluvial2(this)
            %% builds pattern maps for each of the spanning spaces for a cohort, e.g., baseline_cn.

            pwd0 = pushd(this.home);
            for s = asrow(flip(this.selected_spans))
                pwd1 = pushd(fullfile("NumBases"+s, "OPNMF", "niiImg")); 
                this.build_map_for_ggalluvial2(s); 
                popd(pwd1); 
            end
            popd(pwd0);
        end

        function build_map_for_ggalluvial2(this, span)
            %% Identifies time index with rank of basis in the spanning model.

            arguments
                this mladni.NMFHierarchies
                span {mustBeNumeric}
            end
            
            ic = mlfourd.ImagingContext2('Basis_1.nii');
            ic = ic.zeros;
            ifc = ic.nifti;
            for b = 1:span
                ifc1 = mlfourd.ImagingFormatContext2(sprintf('Basis_%i.nii', b));
                ifc.img(:,:,:,b) = ifc1.img;
            end
            
            % identify patterns, arrange to be consistent with table for ggalluvial2
            ifc_argmax = mlfourd.ImagingFormatContext2('Basis_argmax_brain_mask_reordered_patterns.nii.gz');
            ps = unique(ifc_argmax.img);  % pattern ids from 24-model
            ps = ps(ps ~= 0);
            ps = sort(ps);  % ascending
            b = 0;
            for p = asrow(ps)
                b = b + 1;
                hard_cluster = ifc_argmax.img == p;  % max ~ 1
                hard_cluster = reshape(hard_cluster, [1, numel(hard_cluster)]);

                for b1 = 1:span
                    soft_cluster = ifc.img(:,:,:,b1)/max(ifc.img(:,:,:,b1), [], "all");  % max ~ 1
                    soft_cluster = reshape(soft_cluster, [1, numel(soft_cluster)]);
                    similarities(b1) = this.cosine_similarity([hard_cluster; soft_cluster]); %#ok<AGROW>
                end

                [~,selected(b)] = max(similarities); %#ok<AGROW>
            end

            % rearrange ifc
            ifc_final = copy(ifc);
            for b = 1:span
                ifc_final.img(:,:,:,b) = ifc.img(:,:,:,selected(b));
            end
            ifc_final.fileprefix = 'Basis_alluvial'; % composite nifti of bases
            ifc_final.save();
            
            brain_mask = mlfourd.ImagingContext2(fullfile(this.standard_dir, 'MNI152_T1_2mm_brain_mask.nii.gz'));  % tight mask
            ic_final = mlfourd.ImagingContext2(ifc_final);
            ic_final = ic_final .* brain_mask.binarized();
            ic_final.fileprefix = strcat(ifc_final.fileprefix, '_brain_mask');
            ic_final.save();  
        end

        function T = build_table_for_ggalluvial2(this)
            %% maintain voxel identities across pattern-models
            %
            % unique(t.N_Patterns_2)
            %
            % ans =
            %
            %     11
            %     20
            %
            % unique(t.N_Patterns_6)
            %
            % ans =
            %
            %      3
            %     11
            %     13
            %     15
            %     16
            %     20
            %
            % unique(t.N_Patterns_8)
            %
            % ans =
            %
            %      2
            %      3
            %     11
            %     13
            %     15
            %     16
            %     20
            %     24
            %
            % unique(t.N_Patterns_10)
            %
            % ans =
            %
            %      4
            %      5
            %      6
            %      8
            %      9
            %     11
            %     13
            %     16
            %     20
            %     24
            %
            % unique(t.N_Patterns_12)
            %
            % ans =
            %
            %      4
            %      5
            %      6
            %      8
            %      9
            %     11
            %     13
            %     16
            %     18
            %     20
            %     22
            %     24
            %
            % unique(t.N_Patterns_14)
            %
            % ans =
            %
            %      2
            %      3
            %      4
            %      5
            %      6
            %      8
            %      9
            %     11
            %     13
            %     14
            %     16
            %     17
            %     20
            %     24

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
                    img = this.reorder_patterns_by_cosine_similarity(img, img_span_ref);
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

        function heatmap_similarity_vs_yeo7(this)
        end

        function heatmap_similarity_vs_yeo17(this)
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

        function g = pattern_num_from(this, span, basis)
            arguments
                this mladni.NMFHierarchies
                span double
                basis double
            end

            T = this.table_patt_weighted_fdg(span=span);
            g = T{T.indices_bases == basis, 'rank'};
        end

        function img1 = reorder_patterns_by_cosine_similarity(this, img, ref)
            %% Use cosine similarity to reorder pattern indices in img to have best similarity to img_span_ref, e.g. ,
            %  NumBasis24.  Crucial for obtaining interpretable visualization of alluvials.

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

        function T = table_patt_weighted_fdg(this, opts)
            %% reorder bases so that for CN FDG scans, mean(suvr) from component-weighted averaging of sampled scans
            %  is sorted high to low
            %
            % current ans from N=269 :
            %
            % 24×4 table
            %
            %   indices_bases      mu        sigma      rank
            %   _____________    _______    ________    ____
            %
            %        22           1.3443     0.12034      1
            %        10           1.3293     0.11549      2
            %        11           1.3031     0.12068      3
            %        15           1.2762     0.10983      4
            %        13           1.2599     0.10673      5
            %         5           1.2535     0.10701      6
            %        17           1.2461     0.10985      7
            %         9            1.204     0.11401      8
            %         3           1.1831     0.13087      9
            %        16           1.1709    0.094711     10
            %         2           1.1624     0.11145     11
            %        24            1.158    0.097288     12
            %         4           1.1314    0.091663     13
            %        14           1.0844    0.089932     14
            %         8           1.0661    0.071617     15
            %         7           1.0578     0.07942     16
            %        20           1.0444     0.03563     17
            %        12           1.0323     0.12004     18
            %        18           1.0292     0.13535     19
            %         1          0.99408    0.093491     20
            %        23          0.96954     0.06748     21
            %        19          0.93851    0.065269     22
            %        21           0.8152    0.091816     23
            %         6          0.75447    0.048727     24
            %
            % 14×4 table
            %
            %   indices_bases      mu        sigma      rank
            %   _____________    _______    ________    ____
            %
            %         6           1.3253     0.11228      1
            %        10            1.299     0.10874      2
            %         4           1.2347     0.10205      3
            %         3           1.2016       0.124      4
            %        13           1.1986     0.10821      5
            %         2           1.1833    0.097242      6
            %        14           1.1377    0.087597      7
            %         9           1.0418     0.13386      8
            %         7           1.0332     0.11369      9
            %         8           1.0215    0.046537     10
            %         1          0.99491    0.090665     11
            %        11          0.94472     0.06226     12
            %        12          0.87253    0.090973     13
            %         5          0.84685    0.052167     14
            %
            % 12×4 table
            %
            %   indices_bases      mu        sigma      rank
            %   _____________    _______    ________    ____
            %
            %         3           1.3236      0.1146      1
            %        10           1.2893     0.10684      2
            %        11           1.2552     0.10686      3
            %        12           1.1783    0.096855      4
            %         2           1.1772    0.092725      5
            %         7           1.0815      0.1273      6
            %         5           1.0373     0.11823      7
            %         1           0.9917    0.089691      8
            %         6          0.98867    0.041822      9
            %         4          0.95064    0.068452     10
            %         9          0.95016    0.062515     11
            %         8          0.89613    0.086296     12
            %
            % 10×4 table
            %
            %   indices_bases      mu        sigma      rank
            %   _____________    _______    ________    ____
            %
            %         7           1.2936     0.10614      1
            %         9           1.2491     0.10442      2
            %        10           1.2262     0.10252      3
            %         2           1.1429    0.087769      4
            %         5           1.0847     0.12769      5
            %         3           1.0396     0.11667      6
            %         1          0.99766     0.08933      7
            %         4          0.97337    0.041006      8
            %         8          0.95278    0.062537      9
            %         6          0.90491    0.084665     10
            %
            % 8×4 table
            %
            %   indices_bases      mu        sigma      rank
            %   _____________    _______    ________    ____
            %
            %         5           1.2685     0.10098     1
            %         6           1.2294     0.10164     2
            %         4           1.2014     0.11161     3
            %         2           1.0329      0.1192     4
            %         1           1.0299    0.088758     5
            %         3          0.97798    0.042721     6
            %         7          0.96639    0.063746     7
            %         8          0.89032    0.083103     8
            %
            % 6×4 table
            %
            %   indices_bases      mu        sigma      rank
            %   _____________    _______    ________    ____
            %
            %         4           1.2573    0.098568     1
            %         6           1.2052      0.1091     2
            %         2           1.0591     0.11281     3
            %         5           1.0075    0.067739     4
            %         1          0.99994    0.084931     5
            %         3          0.95985    0.041322     6
            %
            % 2×4 table
            %
            %   indices_bases      mu       sigma      rank
            %   _____________    ______    ________    ____
            %
            %         2          1.1019    0.093929     1
            %         1          1.0348    0.069714     2

            arguments
                this mladni.NMFHierarchies
                opts.span = this.N_patterns
            end
            % assert(any(ismember(opts.span, 1:max(this.selected_spans))))  % considering all possible N_patterns

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
        end

        function t = addvars_patterns(this, t_, opts)
            %% Add vars Pattern_1, Pattern_2, to a table that already contains vars Components_1, Components_2, ...
            %  Pattern_* vars are ordered from high to low suvr as prescribed by table_patt_weighted_fdg.

            arguments
                this mladni.NMFHierarchies
                t_ {mustBeNonempty}
                opts.rewrite logical = false
            end

            % perform I/O as needed to load table
            t_fqfp = stackstr();
            if istext(t_) && isfile(t_) && endsWith(t_, ".mat")
                t_fqfp = myfileprefix(t_);
                ld = load(t_);
                t_ = ld.t;
            end
            if istext(t_) && isfile(t_) && endsWith(t_, ".csv'")
                t_fqfp = myfileprefix(t_);
                t_ = readtable(t_, ReadVariableNames=true);
            end
            assert(istable(t_))

            if startsWith(t_.Properties.VariableNames{end}, 'Pattern_')
                t = t_;
                return
            end

            % obtain ordering of components by suvr =: patterns
            tpwf = this.table_patt_weighted_fdg;

            % addvars
            cselect = contains(t_.Properties.VariableNames, "Components");
            c = t_(:, cselect);
            p = c(:, asrow(tpwf.indices_bases));
            p.Properties.VariableNames = strrep(c.Properties.VariableNames, "Components_", "Pattern_");
            t = [t_, p];

            if opts.rewrite
                save(t_fqfp+".mat", "t");
                writetable(t, t_fqfp+".csv", WriteVariableNames=true);
            end
        end
    end

    methods (Static)
        function ARI_between_alluvial_axes(a1, a2)
            pwd0 = pushd(fullfile(getenv("SINGULARITY_HOME"), "ADNI", "NMF_FDG", "baseline_cn"));
            p1 = double(mlfourd.ImagingContext2( ...
                fullfile("NumBases"+a1, "OPNMF", "niiImg", "Basis_argmax_brain_mask.nii")));
            p1 = p1(p1 > 0);
            p2 = double(mlfourd.ImagingContext2( ...
                fullfile("NumBases"+a2, "OPNMF", "niiImg", "Basis_argmax_brain_mask.nii")));
            p2 = p2(p2 > 0);

            ari = mladni.AdjRandIndex.rand_index(p1, p2, "adjusted");
            fprintf("%s: adj. Rand index [%g, %g] -> %g\n", stackstr(), a1, a2, ari);
            popd(pwd0);
        end

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
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
