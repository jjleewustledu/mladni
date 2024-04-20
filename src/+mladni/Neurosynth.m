classdef Neurosynth < handle
    %% Gathers Neurosynth v4-topics-50
    %  
    %  Created 02-Feb-2023 23:20:54 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.13.0.2126072 (R2022b) Update 3 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        N_PATTERNS = mladni.NMF.N_PATTERNS

        % Jones Table 3
        ORDERING = [36 17 12 44 40  41 20 42 22 23  29 47 32 5 27  1 26 25 8 0  19 11 3 48 15  35 49]' + 1
        TERMS = {'Langue comprehension' 'Social' 'Memory' 'Language semantics' 'Negative emotion' ...
            'Visual attention' 'Language perception' 'Numerical' 'Working Memory' 'Emotional cues' ...
            'Reward' 'Response preparation' 'Hearing' 'Facial recognition' 'Addiction' ...
            'Objects' 'Sustenance state' 'Error learning' 'Response inhibition' 'Praxis' ...
            'Stimulus response' 'Motion perception' 'Perception' 'Pain' 'Directed gaze' ...
            'Somatosensory' 'Motor'}'
        PREFIXES27 = {'Langue_comprehension' 'Social' 'Memory' 'Language_semantics' 'Negative_emotion' ...
            'Visual_attention' 'Language_perception' 'Numerical' 'Working_Memory' 'Emotional_cues' ...
            'Reward' 'Response_preparation' 'Hearing' 'Facial_recognition' 'Addiction' ...
            'Objects' 'Sustenance_state' 'Error_learning' 'Response_inhibition' 'Praxis' ...
            'Stimulus_response' 'Motion_perception' 'Perception' 'Pain' 'Directed_gaze' ...
            'Somatosensory' 'Motor'}'
        PREFIXES104 = {'action', 'adaptation', 'addiction', 'anticipation', 'anxiety', 'arousal', 'association', ...
              'autobiographical_memory', 'awareness', 'balance', 'belief', 'categorization', 'choice', ...
              'cognitive_control', 'communication', 'competition', 'concept', 'consciousness', 'consolidation', ...
              'context', 'coordination', 'decision_making', 'discrimination', 'distraction', 'eating', 'efficiency', ...
              'effort', 'emotion_regulation', 'empathy', 'encoding', 'episodic_memory', 'executive_control', ...
              'executive_function', 'expectancy', 'expertise', 'face_recognition', 'facial_expression', ...
              'familiarity', 'fear', 'gaze', 'goal', 'hyperactivity', 'impulsivity', 'induction', 'inference', ...
              'integration', 'intelligence', 'intention', 'interference', 'knowledge', 'language_comprehension', ...
              'learning', 'listening', 'loss', 'maintenance', 'manipulation', 'memory_retrieval', 'mental_imagery', ...
              'monitoring', 'mood', 'morphology', 'motor_control', 'movement', 'multisensory', 'naming', ...
              'navigation', 'object_recognition', 'pain', 'planning', 'priming', 'psychosis', 'reading', ...
              'reasoning', 'recall', 'rehearsal', 'remembering', 'response_inhibition', 'response_selection', ...
              'retention', 'reward', 'rhythm', 'risk', 'rule', 'salience', 'selective_attention', 'semantic_memory', ...
              'sentence_comprehension', 'sleep', 'social_cognition', 'spatial_attention', 'speech_perception', ...
              'speech_production', 'strategy', 'stress', 'sustained_attention', 'thought', 'uncertainty', ...
              'updating', 'valence', 'verbal_fluency', 'visual_attention', 'visual_perception', 'word_recognition', ...
              'working_memory'}
    end

    properties (Dependent)
        bases_volbin
        bases_reported
        imaging
        mask
        workdir
    end

    methods % GET
        function g = get.bases_volbin(this)
            g = 1:24;
        end
        function g = get.bases_reported(this)
            g = this.nmf_radar_.sorted_bases;
        end
        function g = get.imaging(this)
            g = this.imaging_;
        end
        function g = get.mask(this)
            if ~isempty(this.mask_)
                g = this.mask_;
                return
            end
            mask_fqfn = fullfile(getenv('ADNI_HOME'), 'VolBin', 'mask.nii.gz');
            assert(isfile(mask_fqfn))
            this.mask_ = mlfourd.ImagingContext2(mask_fqfn);
            assert(all(size(this.mask_) == size(this.imaging{1})))
            g = this.mask_;
        end
        function     set.mask(this, s)
            this.mask_ = mlfourd.ImagingContext2(s);
        end
        function g = get.workdir(this)
            g = fullfile(getenv('ADNI_HOME'), 'NMF_FDG');
        end
    end
       
    methods
        function T = build_for_brainsmash(this)
            %% Prepares intermediates needed for brainsmash inference described in 
            %  test_neurodegeneration2.py/TestNeurodegeneration2.test_build_stats27.
            %  Requires availability of AFNI.3dmaskdump.
            %
            %  Returns table of 'SummaryTerm', 'TopicTermNumber', 'Filename'}

            workdir = fullfile(getenv('ADNI_HOME'), 'neurosynth.org/v4-topics-50');
            cd(workdir)

            T = table(this.TERMS, this.ORDERING-1);
            T.Properties.VariableNames = {'SummaryTerm', 'TopicTermNumber'};
            Filename = cell(size(T,1), 1);
            AltFileprefix = cell(size(T,1), 1);
            O__ = this.ORDERING;
            T__ = this.TERMS;
            parfor idx = 1:size(T,1)
                idx0 = O__(idx)-1;
                fold = sprintf('topic-%i', idx0);
                fn = glob(fullfile(fold, '*_abs.nii'));
                assert(length(fn) == 1)
                Filename{idx} = fn{1};
                AltFileprefix{idx} = strrep(T__{idx}, ' ', '_');
                system(sprintf('3dmaskdump -mask mask.nii.gz -o %s.txt -noijk %s', AltFileprefix{idx}, Filename{idx}))
            end
            T = addvars(T, Filename, AltFileprefix, NewVariableNames={'Filename', 'AltFileprefix'});
            writetable(T, sprintf('%s_T.csv', stackstr()));
        end
        function T = build_stats_from_logs(this, opts)
            %% e.g.
            % neurodegeneration2_1k =
            % 1648×4 table
            %
            %       r          p           term         basis
            %   _________    _____    ______________    _____
            %
            %     0.03036    0.339    action              1
            %    0.015974     0.64    action              2
            %   -0.083189    0.044    action              3
            %    0.078612    0.026    action              4
            %   -0.043969    0.239    action              5
            %   -0.050962    0.165    action              6
            %   -0.078772    0.033    action              7
            %   -0.063114    0.071    action              8
            %     0.17583        0    action              9
            %    0.099046    0.007    action             10
            %   -0.044881    0.132    action             11
            %    0.081765    0.021    action             12
            %     0.19237        0    action             13
            %       :          :            :             :

            arguments
                this mladni.Neurosynth
                opts.tag {mustBeTextScalar} = "Neurosynth120"
            end

            pwd0 = pushd(fullfile(this.workdir, "brainsmash_output"));

            % call mglob according to opts.tag
            switch convertStringsToChars(opts.tag)
                case 'EB'
                    mg = mglob("EB*.log");
                case 'Neurosynth120'
                    mg = string(this.PREFIXES104) + ".log";
                case 'Neurosynth27'
                    mg = string(this.PREFIXES27) + ".log";
            end
            mg = natsort(mg);

            % parse globbed files, lines from files, then assemble the table
            r = [];
            p = [];
            term = "";
            basis = [];
            for globbed = asrow(mg)
                lines = readlines(globbed(1));
                for a_line = asrow(lines)
                    if isemptytext(a_line(1)) || contains(a_line(1), "Pearson r: nan")
                        continue
                    end
                    try
                        re = regexp(a_line(1), "Pearson r: (?<r>(|-)\d.\d+(|e-\d+))\(p = (?<p>\d.\d+(|e-\d+))\) for (?<term>\S+) basis (?<basis>\d+)", "names");
                        r = [r; str2double(re.r)]; %#ok<*AGROW>
                        p = [p; str2double(re.p)];
                        term = [term; string(re.term)];
                        basis = [basis; str2double(re.basis)];
                    catch 
                        fprintf("File %s failed parsing:  %s\n", globbed(1), a_line(1))
                    end
                end
            end
            term = term(term ~= "");
            T = table(r, p, term, basis);

            save("neurodegeneration2_5k_"+opts.tag+".mat", "T")
            popd(pwd0)
        end
        function [rho,pval] = corr(this, ic, topic, opts)
            arguments
                this mladni.Neurosynth
                ic {mustBeNonempty}
                topic double
                opts.force_zscore = true
            end
            ic = mlfourd.ImagingContext2(ic);
            assert(topic <= length(this.imaging))
            vec1 = ascol(ic.nifti.img(logical(this.mask)));
            if opts.force_zscore
                vec1 = zscore(vec1);
            end
            vec2 = ascol(this.imaging{topic}.nifti.img(logical(this.mask)));
            [rho,pval] = corr([vec1, vec2]);
        end
        function vec = corr_all_topics(this, ic, opts)
            arguments
                this mladni.Neurosynth
                ic {mustBeNonempty}
                opts.force_zscore = true
            end
            vec = nan(50, 1);
            for t = 1:length(this.imaging)
                rho = this.corr(ic, t, force_zscore=opts.force_zscore);
                vec(t) = rho(1,2);
            end
            vec = vec(this.ORDERING);
        end        
        function d = kldiv(this, ic, topic)
            %% See also mlfourd.MatlabTool.kldiv(map, mask)

            arguments
                this mladni.Neurosynth
                ic {mustBeNonempty}
                topic double
            end
            ic = mlfourd.ImagingContext2(ic);
            assert(topic <= length(this.imaging))
            d = ic.kldiv(this.imaging{topic}, this.mask);
        end
        function vec = kldiv_all_topics(this, ic)
            arguments
                this mladni.Neurosynth
                ic {mustBeNonempty}
            end
            vec = nan(50, 1);
            for t = 1:length(this.imaging)
                vec(t) = this.kldiv(ic, t).nifti.img;
            end
            vec = vec(this.ORDERING);
        end
        function s = sim(this, ic, topic, opts)
            %% dot(u,v)/(norm(u)*norm(v))
            %  See also:
            %  https://stackoverflow.com/questions/22432673/how-to-measure-the-cosine-similarity-between-2-images
            %  web(fullfile(docroot, 'stats/pdist.html#mw_39296772-30a1-45f3-a296-653c38875df7'))

            arguments
                this mladni.Neurosynth
                ic {mustBeNonempty}
                topic double
                opts.Distance = 'cosine'
                opts.force_zscore = false
            end
            ic = mlfourd.ImagingContext2(ic);
            assert(topic <= length(this.imaging))
            vec1 = ascol(ic.nifti.img(logical(this.mask)));
            if opts.force_zscore
                vec1 = zscore(vec1);
            end
            vec2 = ascol(this.imaging{topic}.nifti.img(logical(this.mask)));
            %vec1 = decimate(double(vec1), 10);
            %vec2 = decimate(double(vec2), 10);
            s = dot(vec1,vec2)/(norm(vec1)*norm(vec2));
        end
        function vec = sim_all_topics(this, ic, opts)
            arguments
                this mladni.Neurosynth
                ic {mustBeNonempty}
                opts.Distance = 'cosine'
                opts.force_zscore = false
            end
            vec = nan(50, 1);
            for t = 1:length(this.imaging)
                vec(t) = this.sim(ic, t, Distance=opts.Distance, force_zscore=opts.force_zscore);
            end
            vec = vec(this.ORDERING);
        end
        
        function this = Neurosynth(cache, amask)
            arguments
                cache {mustBeText} = ""
                amask = []
            end
            if isfile(cache)
                ld = load(cache);
                this.imaging_ = ld.imaging_cache;
            else
                this.imaging_ = cell(1, 50);
            end
            if ~isempty(amask)
                this.mask_ = mlfourd.ImagingContext2(amask);
            end
            this.j2022_ = mladni.Jones2022();
            this.nmfh_ = mladni.NMFHierarchies();
            this.nmf_radar_ = mladni.NMFRadar();
            this.ns27_ = mladni.Neurosynth27();
            this.ns120_ = mladni.Neurosynth120();

            topics_home = fullfile(getenv('ADNI_HOME'), 'neurosynth.org', 'v4-topics-50');
            pwd0 = pushd(topics_home);            
            for idx = 1:50
                g = glob(sprintf('topic-%i/v4-topics-50_%i_*_association-test_z_FDR_0.01.nii', idx-1, idx-1));
                if ~isempty(g{end})
                    this.imaging_{idx} = mlfourd.ImagingContext2( ...
                        fullfile(topics_home, g{end}));
                end
            end
            imaging_cache = this.imaging_;
            save(fullfile(topics_home, "mladni_Neurosynth_imaging.mat"), "imaging_cache");
            popd(pwd0);
        end        
        
        function T = table_brainsmash_corr(this)
            arguments
                this mladni.Neurosynth
            end

            % FDR correct the corr from external tables
            T_EB = this.j2022_.table;
            T_27 = this.ns27_.table;
            T_120 = this.ns120_.table;
            T_EB.corr(T_EB.fdr > 0.05) = 0;
            T_27.corr(T_27.fdr > 0.05) = 0;
            T_120.corr(T_120.fdr > 0.05) = 0;

            labels_patterns = "P" + (1:24);  % 1:24
            labels_EB = "EB" + (1:10);  % 25:34
            labels_27 = string(asrow(this.ns27_.table.term));  % 35:61
            labels_120 = asrow(this.ns120_.table.term);  % 62:181
            labels = [labels_patterns, labels_EB, labels_27, labels_120];
            N = numel(labels);
            mapping = this.nmfh_.mapping_span_by_suvr;  % reorders patterns to be ranked by their total suvr

            corr = zeros(N, N);
            corr(25:34, 1:24) = T_EB.corr(:, mapping);
            corr(35:61, 1:24) = T_27.corr(:, mapping);
            corr(62:181, 1:24) = T_120.corr(:, mapping);

            T = table(corr, RowNames=asrow(labels));
            T = splitvars(T, 1, NewVariableNames=asrow(labels));

            writetable(T, "Neurosynth_brainsmash_corr.csv", WriteRowNames=true, Delimiter=" ");
        end
        function [T, labels] = table_brainsmash_corr2(this)
            arguments
                this mladni.Neurosynth
            end

            % FDR correct the corr from external tables
            T_27 = this.ns27_.table;
            T_120 = this.ns120_.table;
            T_27.corr(T_27.fdr > 0.05) = 0;
            T_120.corr(T_120.fdr > 0.05) = 0;

            labels_patterns = "P" + (1:24);  % 1:24
            labels_27 = string(asrow(this.ns27_.table.term));  % 25:51
            labels_120 = asrow(this.ns120_.table.term);  % 52:171
            labels = [labels_patterns, labels_27, labels_120];
            N = numel(labels);
            mapping = this.nmfh_.mapping_span_by_suvr;  % reorders patterns to be ranked by their total suvr

            corr = zeros(N, N);
            corr(25:51, 1:24) = T_27.corr(:, mapping);
            corr(52:171, 1:24) = T_120.corr(:, mapping);
            corr = corr';

            T = table(corr, RowNames=asrow(labels));
            T = splitvars(T, 1, NewVariableNames=asrow(labels));

            writetable(T, "Neurosynth_brainsmash_corr2.csv", WriteRowNames=true, Delimiter=" ");
        end
        function T = table3(this)

            T = table(this.TERMS, this.ORDERING);
            T.Properties.VariableNames = {'Summary term', 'Topic term number'};

            pwd0 = pushd(sprintf('%s/NMF_FDG/baseline_cn/NumBases%i/OPNMF/niiImg', getenv('ADNI_HOME'), this.N_PATTERNS));
            for idx = 1:this.N_PATTERNS
                vec = this.corr_all_topics(sprintf('Basis_%i.nii', idx));
                T = addvars(T, vec);
                T.Properties.VariableNames{end} = sprintf('P%i', idx);
            end

            cd('/Volumes/PrecunealSSD/Singularity/ADNI/neurovault.org/gradient_1_cortical_')
            % this.mask = 'volume.cort.0.nii.gz';
            FC_gradient = this.corr_all_topics('volume.cort.0.nii.gz');
            T = addvars(T, FC_gradient);
            T.Properties.VariableNames{end} = 'FC gradient';
            popd(pwd0)
        end
        function h = heatmap3(this)
            %% Pearson correlations

            pwd0 = pushd(sprintf('%s/NMF_FDG/baseline_cn/NumBases%i/OPNMF/niiImg', getenv('ADNI_HOME'), this.N_PATTERNS));
            mat = [];
            for idx = 1:this.N_PATTERNS
                mat = [mat, this.corr_all_topics(sprintf('Basis_%i.nii', idx))]; %#ok<AGROW> 
            end
            cd('/Volumes/PrecunealSSD/Singularity/ADNI/neurovault.org/gradient_1_cortical_')
            % this.mask = 'volume.cort.0.nii.gz';
            mat = [mat, this.corr_all_topics('volume.cort.0.nii.gz')];

            clbls = {'P1' 'P2' 'P3' 'P4' 'P5' 'P6' 'P7' 'P8' 'P9' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'FC Gradient'};
            rlbls = this.TERMS;
            h = mladni.Jones2022.heatmap(mat, clbls, rlbls);
            ylabel('Summary Terms')
            xlabel('ADNI Pattern, FC Gradient')

            popd(pwd0);
        end
        function h = heatmap3_cossim(this)
            %% Cosine similarity

            pwd0 = pushd('/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_cn/NumBases16/OPNMF/niiImg');
            mat = [];
            for idx = 1:this.N_PATTERNS
                mat = [mat, this.sim_all_topics(sprintf('Basis_%i.nii', idx))]; %#ok<AGROW> 
            end
            cd('/Volumes/PrecunealSSD/Singularity/ADNI/neurovault.org/gradient_1_cortical_')
            this.mask = 'volume.cort.0.nii.gz';
            mat = [mat, this.sim_all_topics('volume.cort.0.nii.gz')];

            clbls = {'P1' 'P2' 'P3' 'P4' 'P5' 'P6' 'P7' 'P8' 'P9' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'FC Gradient'};
            rlbls = this.TERMS;
            h = mladni.Jones2022.heatmap(mat, clbls, rlbls);
            ylabel('Summary Terms')
            xlabel('ADNI Pattern, FC Gradient')

            popd(pwd0);
        end
        function h = heatmap3_kldiv(this)
            %% Kuehback-Leibler divergence, relative entropy from topical-term maps to patterns

            pwd0 = pushd('/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_cn/NumBases16/OPNMF/niiImg');
            mat = [];
            for idx = 1:this.N_PATTERNS
                mat = [mat, this.kldiv_all_topics(sprintf('Basis_%i.nii', idx))]; %#ok<AGROW> 
            end
            cd('/Volumes/PrecunealSSD/Singularity/ADNI/neurovault.org/gradient_1_cortical_')
            this.mask = 'volume.cort.0.nii.gz';
            mat = [mat, this.kldiv_all_topics('volume.cort.0.nii.gz')];

            clbls = {'P1' 'P2' 'P3' 'P4' 'P5' 'P6' 'P7' 'P8' 'P9' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'FC Gradient'};
            rlbls = this.TERMS;
            h = mladni.Jones2022.heatmap(mat, clbls, rlbls, CellLabelFormat='%3.0f', ColorScaling='scaledrows', FlipColormap=true);
            ylabel('Summary Terms')
            xlabel('ADNI Pattern, FC Gradient')

            popd(pwd0);
        end
    end

    methods (Static)
        function add_all()
            ga = glob(fullfile('topic-*', 'v4-topics-50_*association*.nii'))';
            gu = glob(fullfile('topic-*', 'v4-topics-50_*uniformity*.nii'))';
            assert(length(ga) == length(gu))
            N = length(ga);

            assoc = mlfourd.ImagingContext2(ga{1});
            unif = mlfourd.ImagingContext2(gu{1});
            for n = 2:N
                assoc = assoc + mlfourd.ImagingContext2(ga{n});
                assoc.fileprefix = 'v4-topics-50_association-test_z_FDR_0.01';
                unif = unif + mlfourd.ImagingContext2(gu{n});
                unif.fileprefix = 'v4-topics-50_uniformity-test_z_FDR_0.01';
            end
            assoc.save
            unif.save
        end
        function view_qc(varargin)
            added = mlfourd.ImagingContext2( ...
                fullfile(getenv('ADNI_HOME'), 'neurosynth.org', 'v4-topics-50', ...
                'v4-topics-50_association-test_z_FDR_0.01.nii')); % product of this.add_all()
            added.view_qc(varargin{:})
        end
    end

    %% PRIVATE

    properties (Access = private)
        imaging_
        j2022_
        mask_
        nmfh_
        nmf_radar_
        ns27_
        ns120_
    end

    methods (Access = private)
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
