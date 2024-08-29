classdef Topology < handle
    %% line1
    %  line2
    %  
    %  Created 29-Jul-2024 22:50:15 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 24.1.0.2653294 (R2024a) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    methods
        function ic = mask(this)
            if ~isempty(this.mask_)
                ic = this.mask_;
                return
            end

            ic = mlfourd.ImagingContext2( ...
                fullfile(getenv("SINGULARITY_HOME"), "ADNI", "VolBin", "mask.nii.gz"));
            this.mask_ = ic;
        end
        function d = niiImgDir(this, opts)
            arguments
                this mladni.Topology
                opts.num_bases int16 = 24
            end

            d = fullfile(getenv("SINGULARITY_HOME"), ...
                "ADNI", "NMF_FDG", "baseline_cn", "NumBases" + opts.num_bases, "OPNMF", "niiImg");
        end
        function ic = patterns2logp(this, opts)
            %% Open sets in the topological space will be represented by prob over the MNI brain,
            %  renornalized by uniform distrib. ~ 1/numel(mask)
            %  Pattern indices are by rank of SUVR of pattern-weighted averages of ADNI FDG.

            arguments
                this mladni.Topology
                opts.num_bases int16 = 24
                opts.do_checks logical = true
                opts.use_cache logical = true
            end

            if opts.use_cache
                fqfn = fullfile(this.niiImgDir(num_bases=opts.num_bases), stackstr() + ".nii");
                if isfile(fqfn)
                    ic = mlfourd.ImagingContext2(fqfn);
                    return
                end

                % no filesystem cache available, so build from Basis_alluvial_brain_mask.nii
            end

            all = mlfourd.ImagingContext2( ...
                fullfile(this.niiImgDir(num_bases=opts.num_bases), "Basis_alluvial_brain_mask.nii"));
            % all = all.masked(this.mask);

            all_ifc = all.imagingFormat;
            for t = 1:opts.num_bases
                all_ifc.img(:,:,:,t) = log(all_ifc.img(:,:,:,t)) - log(sum(all_ifc.img(:,:,:,t), "all"));
            end
            ic = mlfourd.ImagingContext2(all_ifc);
            ic.fileprefix = stackstr(use_underscores=true);
            ic.save();

            if opts.do_checks
                this.check_sum_max(all_ifc);
            end
        end
        function ic = patterns2probabilities(this, opts)
            %% Open sets in the topological space will be represented by prob over the MNI brain;
            %  Pattern indices are by rank of SUVR of pattern-weighted averages of ADNI FDG.

            arguments
                this mladni.Topology
                opts.num_bases int16 = 24
                opts.do_checks logical = true
                opts.use_cache logical = true
            end

            if opts.use_cache
                fqfn = fullfile(this.niiImgDir(num_bases=opts.num_bases), stackstr() + ".nii");
                if isfile(fqfn)
                    ic = mlfourd.ImagingContext2(fqfn);
                    return
                end
                if ~isempty(this.pattern_probabilities_)
                    ic = copy(this.pattern_probabilities_);
                    return
                end

                % no cache available, so build this.pattern_probabilities_
            end

            all = mlfourd.ImagingContext2( ...
                fullfile(this.niiImgDir(num_bases=opts.num_bases), "Basis_all.nii"));
            all = all.masked(this.mask);

            all_ifc = all.imagingFormat;
            for t = 1:opts.num_bases
                all_ifc.img(:,:,:,t) = all_ifc.img(:,:,:,t)/sum(all_ifc.img(:,:,:,t), "all");
            end
            ic = mlfourd.ImagingContext2(all_ifc);
            ic.fileprefix = stackstr(use_underscores=true);
            this.pattern_probabilities_ = ic;

            if opts.do_checks
                this.check_sum_max(all_ifc);
            end
        end
        function ic = patterns2prob_renorm(this, opts)
            %% Open sets in the topological space will be represented by prob over the MNI brain,
            %  renornalized by uniform distrib. ~ 1/numel(mask)
            %  Pattern indices are by rank of SUVR of pattern-weighted averages of ADNI FDG.

            arguments
                this mladni.Topology
                opts.num_bases int16 = 24
                opts.do_checks logical = true
                opts.use_cache logical = true
            end

            if opts.use_cache
                fqfn = fullfile(this.niiImgDir(num_bases=opts.num_bases), stackstr() + ".nii");
                if isfile(fqfn)
                    ic = mlfourd.ImagingContext2(fqfn);
                    return
                end

                % no filesystem cache available, so build from Basis_alluvial_brain_mask.nii
            end

            all = mlfourd.ImagingContext2( ...
                fullfile(this.niiImgDir(num_bases=opts.num_bases), "Basis_alluvial_brain_mask.nii"));
            % all = all.masked(this.mask);

            all_ifc = all.imagingFormat;
            Nmask = dipsum(this.mask.binarized);  % ~228483
            for t = 1:opts.num_bases
                all_ifc.img(:,:,:,t) = Nmask*all_ifc.img(:,:,:,t)/sum(all_ifc.img(:,:,:,t), "all");
            end
            ic = mlfourd.ImagingContext2(all_ifc);
            ic.fileprefix = stackstr(use_underscores=true);
            ic.save();

            if opts.do_checks
                this.check_sum_max(all_ifc);
            end
        end
        function T = table_overlap_with_Yeo7(this, obj, yeo_index)
            arguments
                this mladni.Topology
                obj {mustBeNonempty} = "Topology_patterns2probabilities_1mm.nii.gz"
                yeo_index double = 7  % DMN
            end

            ic = mlfourd.ImagingContext2(obj);
            if ~contains(ic.fileprefix, "_1mm")
                ic = mlfourd.ImagingContext2(this.flirt(ic.fqfn));
                assert(contains(ic.fileprefix, "_1mm"))
            end
            assert(isfile(ic.fqfn))

            msk = mlfourd.ImagingContext2( ...
                fullfile(getenv("REFDIR"), "MNI152_T1_1mm_brain_mask.nii.gz"));
            yeodir = fullfile(getenv("REFDIR"), "Yeo");
            yeo7 = mlfourd.ImagingContext2( ...
                fullfile(yeodir, "Yeo2011_7Networks_to_MNI152_T1_1mm.nii.gz"));

            yeo_img = single(yeo7.imagingFormat.img == yeo_index);
            yeo_img = asrow(yeo_img(:));
            ic = ic.masked(msk);
            pattern = ascol(1:size(ic, 4));
            similarity = nan(size(ic, 4), 1);
            for t = asrow(pattern)
                obj_img = ic.imagingFormat.img(:,:,:,t);
                obj_img = asrow(obj_img(:));
                similarity(t) = this.cosine_similarity([yeo_img; obj_img]);
            end

            T = table(pattern, similarity);
            T = sortrows(T, "similarity", "descend");
        end
        function ic = union_patterns(this, opts)
            %% nerve is homotopic to union of open subsets (Ghrist, Elem. Appl. Top., pg. 31)

            arguments
                this mladni.Topology
                opts.except double = []
            end

            ic = this.patterns2probabilities();
            if ~isempty(opts.except)
                ifc = ic.imagingFormat;
                ifc.img(:,:,:, opts.except) = [];
                ic = mlfourd.ImagingContext2(ifc);
            end
            ic = timeAveraged(ic); 
            ic = ic/sum(ic, "all");
            ic.fileprefix = stackstr(use_underscores=true);            
            if ~isempty(opts.except)
                ic.fileprefix = ic.fileprefix + "_except" + opts.except;
            end
        end
        function view_nerves(this, opts)
            arguments
                this mladni.Topology
                opts.save_niis logical = true
            end

            model_degrees = [2, 6, 8, 10, 12, 14, 24];

            md_idx = 0;
            for model_degree = model_degrees(1:end-1)
                md_idx = md_idx + 1;

                nodes = this.patterns2logp(num_bases=model_degree, do_checks=true);
                % fprintf("%s:  view nodes of %i-model\n", stackstr(), model_degree);
                % nodes.view()

                nodes_next = this.patterns2logp(num_bases=model_degrees(md_idx+1), do_checks=true);
                % fprintf("%s:  view nodes of %i-model\n", stackstr(), model_degrees(md_idx+1));
                % nodes_next.view()
                
                for node_index = 1:model_degree
                    edges = this.node2edges(nodes, node_index=node_index, nodes_next=nodes_next);
                    fprintf("%s:  view %i edges of node %i of %i-model\n", stackstr(), size(edges, 4), node_index, model_degree);
                    % edges.view()
                    if opts.save_niis
                        edges.save(); 
                    end
                end
                
                edges_sumt = mlfourd.ImagingContext2(edges);
                edges_sumt = edges_sumt.timeAveraged();
                edges_sumt.view();
                if opts.save_niis
                    edges_sumt.save();
                end
            end
        end
        function sim_1 = view_dominant_paths(this)
            %% t = mladni.Topology; t.view_dominant_paths
            %
            % t =
            %
            %   Topology with no properties.
            %
            % model degree: 2, alluvial id 11, best next-alluvial id 11, best sim 0.903965
            % model degree: 2, alluvial id 20, best next-alluvial id 20, best sim 0.898907
            %
            % model degree: 6, alluvial id 3, best next-alluvial id 3, best sim 0.989825
            % model degree: 6, alluvial id 11, best next-alluvial id 11, best sim 0.974347
            % model degree: 6, alluvial id 13, best next-alluvial id 13, best sim 0.98936
            % model degree: 6, alluvial id 15, best next-alluvial id 15, best sim 0.987104
            % model degree: 6, alluvial id 16, best next-alluvial id 16, best sim 0.984398
            % model degree: 6, alluvial id 20, best next-alluvial id 20, best sim 0.924761
            %
            % model degree: 8, alluvial id 2, best next-alluvial id 5, best sim 0.992126  <=
            % model degree: 8, alluvial id 3, best next-alluvial id 8, best sim 0.972389  <=
            % model degree: 8, alluvial id 11, best next-alluvial id 11, best sim 0.972525
            % model degree: 8, alluvial id 13, best next-alluvial id 13, best sim 0.993593
            % model degree: 8, alluvial id 15, best next-alluvial id 6, best sim 0.990697  <=
            % model degree: 8, alluvial id 16, best next-alluvial id 16, best sim 0.980674
            % model degree: 8, alluvial id 20, best next-alluvial id 20, best sim 0.985872
            % model degree: 8, alluvial id 24, best next-alluvial id 24, best sim 0.990014
            %
            % model degree: 10, alluvial id 4, best next-alluvial id 18, best sim 0.990785  <=
            % model degree: 10, alluvial id 5, best next-alluvial id 5, best sim 0.995395
            % model degree: 10, alluvial id 6, best next-alluvial id 6, best sim 0.991111
            % model degree: 10, alluvial id 8, best next-alluvial id 8, best sim 0.984516
            % model degree: 10, alluvial id 9, best next-alluvial id 9, best sim 0.995649
            % model degree: 10, alluvial id 11, best next-alluvial id 11, best sim 0.994069
            % model degree: 10, alluvial id 13, best next-alluvial id 13, best sim 0.995966
            % model degree: 10, alluvial id 16, best next-alluvial id 16, best sim 0.996291
            % model degree: 10, alluvial id 20, best next-alluvial id 20, best sim 0.994727
            % model degree: 10, alluvial id 24, best next-alluvial id 24, best sim 0.993142
            %
            % model degree: 12, alluvial id 4, best next-alluvial id 14, best sim 0.959376  <=
            % model degree: 12, alluvial id 5, best next-alluvial id 2, best sim 0.987929  <=
            % model degree: 12, alluvial id 6, best next-alluvial id 6, best sim 0.983694
            % model degree: 12, alluvial id 8, best next-alluvial id 8, best sim 0.991612
            % model degree: 12, alluvial id 9, best next-alluvial id 9, best sim 0.97357
            % model degree: 12, alluvial id 11, best next-alluvial id 11, best sim 0.985669
            % model degree: 12, alluvial id 13, best next-alluvial id 13, best sim 0.99485
            % model degree: 12, alluvial id 16, best next-alluvial id 16, best sim 0.99596
            % model degree: 12, alluvial id 18, best next-alluvial id 4, best sim 0.993759  <=
            % model degree: 12, alluvial id 20, best next-alluvial id 20, best sim 0.991975
            % model degree: 12, alluvial id 22, best next-alluvial id 17, best sim 0.992946  <=
            % model degree: 12, alluvial id 24, best next-alluvial id 24, best sim 0.978966
            %
            % model degree: 14, alluvial id 2, best next-alluvial id 2, best sim 0.974313
            % model degree: 14, alluvial id 3, best next-alluvial id 3, best sim 0.976274
            % model degree: 14, alluvial id 4, best next-alluvial id 18, best sim 0.9694  <=
            % model degree: 14, alluvial id 5, best next-alluvial id 5, best sim 0.972536 
            % model degree: 14, alluvial id 6, best next-alluvial id 15, best sim 0.962829  <=
            % model degree: 14, alluvial id 8, best next-alluvial id 8, best sim 0.974189
            % model degree: 14, alluvial id 9, best next-alluvial id 9, best sim 0.97552
            % model degree: 14, alluvial id 11, best next-alluvial id 11, best sim 0.972337
            % model degree: 14, alluvial id 13, best next-alluvial id 13, best sim 0.97208
            % model degree: 14, alluvial id 14, best next-alluvial id 14, best sim 0.972822
            % model degree: 14, alluvial id 16, best next-alluvial id 16, best sim 0.972373
            % model degree: 14, alluvial id 17, best next-alluvial id 17, best sim 0.966144
            % model degree: 14, alluvial id 20, best next-alluvial id 20, best sim 0.958671
            % model degree: 14, alluvial id 24, best next-alluvial id 24, best sim 0.968645
            %
            % ans =
            %
            %   6×1 cell array
            %
            %     { 2×6  double}
            %     { 6×8  double}
            %     { 8×10 double}
            %     {10×12 double}
            %     {12×14 double}
            %     {14×24 double}

            bin_mask = logical(this.mask);
            model_degrees = [2, 6, 8, 10, 12, 14, 24];

            sim_1 = cell(numel(model_degrees) - 1, 1);  % cells of models ~ {m_2, m_6, ..., m_24}
            md_idx = 0;
            for model_degree = model_degrees(1:end-1)
                md_idx = md_idx + 1;

                nodes = this.patterns2logp(num_bases=model_degree, do_checks=true);
                % fprintf("%s:  view nodes of %i-model\n", stackstr(), model_degree);
                % nodes.view()

                nodes_next = this.patterns2logp(num_bases=model_degrees(md_idx+1), do_checks=true);
                % fprintf("%s:  view nodes of %i-model\n", stackstr(), model_degrees(md_idx+1));
                % nodes_next.view()
                
                % use interpretable alluvial IDs p1, ..., p24, with "hue" colormapping
                aids = this.model_degree_to_alluvial_ids(model_degree);
                next_aids = this.model_degree_to_alluvial_ids(model_degrees(md_idx+1));

                sim_2 = nan(model_degree, model_degrees(md_idx+1));  % matrix of similarities:  model_m -> mode_{m+1}
                for node_index = 1:model_degree
                    node = nodes.imagingFormat.img(:,:,:,node_index);
                    node = asrow(node(bin_mask));
                    
                    for node_next_index = 1:model_degrees(md_idx+1)
                        node_next = nodes_next.imagingFormat.img(:,:,:,node_next_index);
                        node_next = asrow(node_next(bin_mask));
                        
                        sim_2(node_index, node_next_index) = this.cosine_similarity([node; node_next]);
                    end
                    [best_sim, best_idx] = max(sim_2(node_index, :), [], 2);
                    fprintf("model degree: %g, alluvial id %i, best next-alluvial id %i, best sim %g\n", ...
                        model_degree, aids(node_index), next_aids(best_idx), best_sim);
                end
                
                sim_1{md_idx} = sim_2;
            end

            %T = table();
        end
        function ic = all_intersections_patterns(this)
            %% a sheaf subordinate to the cover of open sets assigns an abelian group to every non-empty 
            %  intersection of cover elements (Ghrist, Elem. Appl. Top., pg. 190);
            %  Loring Tu on graph of a good cover:  https://www.youtube.com/watch?v=BK3VO_1jWqA @ 2:30
            %  "de Rham cohomology depens only on the simplicial complex generated by the Nerve"

            pp = this.patterns2probabilities();
            pp_ifc = pp.imagingFormat;
            Nt = size(pp, 4);

            sz = size(pp);
            img_accum = zeros(sz(1:3));
            for t = 1:Nt
                img_t = pp_ifc.img(:,:,:,t);  % normalized
                others = this.union_patterns(except=t);
                img_others = others.imagingFormat.img;  % normalized
                img_accum = img_accum + img_t .* img_others;
            end
            img_accum = img_accum/sum(img_accum, "all");

            ifc = copy(pp_ifc);
            ifc.img = img_accum;
            ic = mlfourd.ImagingContext2(ifc);
            ic.fileprefix = stackstr(use_underscores=true);
        end
    end

    %% PRIVATE

    properties (Access = private)
        mask_
        pattern_probabilities_
    end

    methods (Access = private)
        function check_sum_max(~, all_ifc)
            fprintf("%s: checking pattern sums\n", stackstr());
            for t = 1:size(all_ifc, 4)
                img = all_ifc.img(:,:,:,t);
                img = img(isfinite(img));
                fprintf("sum_%i:  %g\n", t, sum(img, "all")); 
            end

            fprintf("%s: checking pattern maxes\n", stackstr());
            for t = 1:size(all_ifc, 4)
                img = all_ifc.img(:,:,:,t);
                img = img(isfinite(img));
                fprintf("max_%i:  %g\n", t, max(img, [], "all")); 
            end
        end
        function sim = cosine_similarity(~, A)
            sim = 1 - pdist(A, 'cosine');
        end
        function out = flirt(~, in)
            flirt = fullfile(getenv("FSLDIR"), "bin", "flirt");
            out = myfileprefix(in) + "_1mm.nii.gz";
            init = "/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_cn/NumBases24/OPNMF/niiImg/MNI152_T1_2mm_on_1mm.mat";
            ref = "/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_cn/NumBases24/OPNMF/niiImg/MNI152_T1_1mm.nii.gz";
            cmd = sprintf("%s -in %s -applyxfm -init %s -out %s -paddingsize 0.0 -interp spline -ref %s", ...
                flirt, in, init, out, ref);
            system(cmd)
        end
        function edges = node2edges(this, nodes, opts)
            %% returns edges:  log P(nodes_next, nodes(node_index))

            arguments
                this mladni.Topology
                nodes mlfourd.ImagingContext2 {mustBeNonempty}
                opts.node_index double = 1
                opts.nodes_next mlfourd.ImagingContext2 {mustBeNonempty}
            end
            
            node = nodes.imagingFormat;
            node.img = node.img(:,:,:, opts.node_index);

            nodes_next = opts.nodes_next.imagingFormat;

            edges = copy(nodes_next);
            Nedges = size(nodes_next, 4);
            for e = 1:Nedges
                edges.img(:,:,:,e) = node.img + nodes_next.img(:,:,:,e);
            end
            edges.fileprefix = nodes.fileprefix + "_node" + opts.node_index + "to" + Nedges + "edges";
        end
        function ids = model_degree_to_alluvial_ids(~, model_degree)
            switch model_degree
                case 2
                    ids = [11,20];
                case 6
                    ids = [3,11,13,15,16,20];
                case 8
                    ids = [2,3,11,13,15,16,20,24];
                case 10
                    ids = [4,5,6,8,9,11,13,16,20,24];
                case 12
                    ids = [4,5,6,8,9,11,13,16,18,20,22,24];
                case 14
                    ids = [2,3,4,5,6,8,9,11,13,14,16,17,20,24];
                case 24
                    ids = 1:24;
                otherwise
                    error("mladni:ValueError", stackstr())
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
