classdef Neurosynth120 < handle
    %% Builds analog of Alexander-Bloch et al Neuroimage 178 (2018) 540.
    %  https://doi.org/10.1016/j.neuroimage.2018.05.070.
    %  
    %  Created 26-Apr-2023 23:05:11 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.14.0.2239454 (R2023a) Update 1 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        nbases
        neurosynthdir
        workdir
    end

    properties (Dependent)
    end

    methods
        function this = Neurosynth120(nbases)
            arguments
                nbases double {mustBeScalarOrEmpty} = 16
            end
            this.nbases = nbases;
            this.neurosynthdir = fullfile(getenv('ADNI_HOME'), 'neurosynth.org');
            this.workdir = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
        end

        function disp_even_labels(this)
            T = this.table();
            disp(T.term(mod(T.index1,2) == 0))
            disp(T.groups(mod(T.index1,2) == 0))
        end
        function disp_odd_labels(this)
            T = this.table();
            disp(T.term(mod(T.index1,2) == 1))
            disp(T.groups(mod(T.index1,2) == 1))
        end
        function h = heatmap(this, mat, clbls, rlbls, opts)
            arguments
                this mladni.Neurosynth120
                mat double = []
                clbls cell = {}
                rlbls cell = {}
                opts.CellLabelFormat {mustBeTextScalar} = '%0.2g'
                opts.Colormap double = viridis
                opts.ColorScaling {mustBeTextScalar} = 'scaled' % 'scaledrows' 'scaledcolumns'
                opts.FlipColormap logical = false
                opts.FontSize {mustBeScalarOrEmpty} = 18
            end
            if isempty(mat)
                T = this.table();
                mat = T.corr;
                mat(T.fdr > 0.05) = 0;
            end
            if isempty(clbls)
                clbls = num2cell(1:this.nbases);
                clbls = cellfun(@(x) sprintf('ADNI P%i', x), clbls, UniformOutput=false);
            end
            if isempty(rlbls)
                rlbls = asrow(convertStringsToChars(T.term));
            end
            if opts.FlipColormap
                opts.Colormap = flipud(opts.Colormap);
            end
            figure;
            h = heatmap(clbls, rlbls, mat, ...
                CellLabelColor='none', ...
                CellLabelFormat=opts.CellLabelFormat, ...
                Colormap=opts.Colormap, ColorScaling=opts.ColorScaling, ...
                FontSize=opts.FontSize);
            h.Title = "Pearson correlation of 104 Neurosynth terms and ADNI patterns";
        end
        function patterns_for_term(this, term)
            T = this.table();
            U = T(T.term == term, :);
            found = find(U.fdr <= 0.05);
            if ~isempty(found)
                fprintf("correlation with %s significant for pattern: %s\n", term, mat2str(found))
            end
        end
        function patterns_for_terms(this)
            terms = {'autobiographical memory', 'semantic memory', 'recall', 'memory retrieval', 'encoding', 'remembering', 'episodic memory', ...
                'coordination', 'motor control', 'movement', 'planning', 'action', 'representation', ...
                'choice', 'decision making', 'reward', 'addiction', 'anticipation', ...
                'arousal', 'valence', 'fear', 'anxiety', ...
                'search', 'visual search', ...
                'interference', 'working memory', 'maintenance', ...
                'spatial attention', 'selective attention', ...
                'speech perception', 'speech production', 'listening', 'language comprehension', 'sentence comprehension', 'meaning', 'reading', 'naming', 'word recognition'};
            for t = terms
                this.patterns_for_term(t{1});
            end
        end
        function T = table_built_stats(this, varargin)
            %% test_neurodegeneration2.test_build_stats;
            %  hand assembled by importing neurodegeneration2_1k.log

            if ~isempty(this.table_built_stats_)
                T = this.table_built_stats_;
                return
            end
            ld = load(fullfile(this.neurosynthdir, 'neurodegeneration2_1k.mat'));
            this.table_built_stats_ = ld.neurodegeneration2_1k;
            T = this.table_built_stats_;            
            T = this.table_paren(T, varargin{:});
        end
        function T = table_termlist(this, varargin)
            %% enumerations of Alexander-Bloch Fig. 2a
            %  hand assembled from termlist.numbers

            if ~isempty(this.table_termlist_)
                T = this.table_termlist_;
                return
            end
            ld = load(fullfile(this.neurosynthdir, 'termlist.mat'));
            this.table_termlist_ = ld.termlist;
            T = this.table_termlist_;            
            T = this.table_paren(T, varargin{:});
        end
        function T = table(this, varargin)
            %% compendium
            %  VariableNames:  'index', 'odd', 'term', 'corr', 'pval', 'fdr', 'groups'

            if ~isempty(this.table_)
                T = this.table_;
                return
            end

            % aufbau new table variables
            this.table_ = this.table_termlist;
            corr = zeros(length(this.table_.index), this.nbases);
            pval = ones(length(this.table_.index), this.nbases);
            fdr = ones(length(this.table_.index), this.nbases);
            this.table_ = addvars(this.table_, corr, pval, fdr, 'After', 'term', ...
                'NewVariableNames', {'corr', 'pval', 'fdr'});

            terms = unique(this.table_built_stats.term); % categorical
            for it = 1:length(terms)
                term_ = terms(it);

                % find idx of term_ in this.table_termlist
                it1 = find(this.table_.term == string(term_)); 
                if isempty(it1)
                    continue
                end
                assert(isscalar(it1))

                % term_-selected subtable of this.table_built_stats; 
                % avoid caches
                S = this.table_built_stats_(this.table_built_stats_.term == term_, ':'); 
                S = sortrows(S, 'basis');
                if this.nbases == length(S.basis)
                    r_row = S.r';
                    p_row = S.p';
                    this.table_.corr(it1,:) = r_row;
                    this.table_.pval(it1,:) = p_row;
                    continue
                end

                % special cases of U.basis
                r_row = nan(1, this.nbases);
                p_row = nan(1, this.nbases);
                for abasis = S.basis % int-valued
                    r_row(abasis) = S.r(abasis);
                    p_row(abasis) = S.p(abasis);
                end
                this.table_.corr(it1,:) = r_row;
                this.table_.pval(it1,:) = p_row;
            end
            
            % trim new table; update index -> index1
            %this.table_ = this.table_(~contains(this.table_.filename, ".txt"), ':');
            %this.table_ = this.table_(~any(isnan(this.table_.corr), 2), ':');
            %this.table_ = this.table_(~any(isnan(this.table_.pval), 2), ':');
            Nrows = size(this.table_, 1);
            this.table_ = addvars(this.table_, (1:Nrows)', 'After', 'index', 'NewVariableNames', 'index1');

            % do fdr
            for bi = 1:this.nbases
                [~,~,~,P] = fdr_bh(this.table_.pval(:,bi), 0.05, 'dep', 'yes');
                this.table_.fdr(:,bi) = ascol(P);
            end

            T = this.table_;
            T = this.table_paren(T, varargin{:});
        end
    end

    methods (Static)
        function t = table_paren(varargin)
            t = mladni.AdniMerge.table_paren(varargin{:});
        end
    end

    %% PRIVATE

    properties (Access = private)
        table_
        table_built_stats_
        table_termlist_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
