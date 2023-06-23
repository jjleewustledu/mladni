classdef Anticlust < handle
    %% ANTICLUST
    %     R snippets based on Petra's snippets:
    %     =====================================
    %     install.packages("anticlust")
    %     library(anticlust)
    %     install.packages("dplyr")
    %     library(dplyr)
    %     setwd("/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_pve1_msk/NumBases16/components")
    %     subjects = read.csv("anticlust.csv")
    %     continuous.vars <- subjects[, c("MergeAge")]
    %     categorical.vars <- subjects[, c("MergePtGender", "AmyloidStatus", "CDGLOBAL")]
    %     v <- anticlustering(continuous.vars, 2, categories=categorical.vars)
    %     subjects$Split <- v
    %     repA = subjects[subjects$Split == 1, ]
    %     repB = subjects[subjects$Split == 2, ]
    %     write.table(repA, "anticlust_repA.csv", row.names = F, col.names = T, sep = ',', quote = F)
    %     write.table(repB, "anticlust_repB.csv", row.names = F, col.names = T, sep = ',', quote = F)
    %  
    %  Created 18-Aug-2022 21:27:03 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2009381 (R2022a) Update 4 for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        label
        repA
        repB
        workpath
    end

    properties (Constant)
        subgroups0 = { ...
                'cn', 'preclinical', ...
                'cdr_0p5_apos_emci', 'cdr_0p5_apos_lmci', 'cdr_0p5_apos_mci', 'cdr_gt_0p5_apos'}
        subgroups = { ...
            'cn', 'preclinical', 'cdr_0p5_apos', 'cdr_gt_0p5_apos', 'cdr_ge_0p5_aneg'}

    end

    properties (Dependent)
        folder_repA
        folder_repB
    end

    methods % GET
        function g = get.folder_repA(this)
            g = sprintf('baseline_%s_repA', this.label);
        end
        function g = get.folder_repB(this)
            g = sprintf('baseline_%s_repB', this.label);
        end
    end

    methods (Static)
        function T = addvars_Filename(T0)
            if ~istable(T0) && isfile(T0)
                if contains(T0, '.csv')
                    T0 = readtable(T0);
                else
                    ld = load(T0);
                    T0 = ld.(mybasename(T0));
                end
            end
            assert(any(contains(T0.Properties.VariableNames, 'Subject')))

            Filename = [];
            for r = 1:size(T0, 1)
                try
                    re = regexp(T0.Subject{r}, '(?<pre>\d{3})_S_(?<post>\d{4})', 'names');
                    sub = sprintf('sub-%sS%s', re.pre, re.post);
                    globbing = glob(fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', sub, 'ses-*')); 
                    ses = basename(globbing{1});
                    globbing = glob(fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives', sub, ses, 'pet', 'sub-*S*_ses-*_trc-FDG_proc-CASU*_orient-rpi_pet_on_T1w_Warped_dlicv.nii.gz'));
                    % with or without ponsvermis renorm
                    Filename = [Filename; string(globbing{1})]; %#ok<AGROW> 
                catch 
                    fprintf('defective files for %s\n', T0.Subject{r})
                    Filename = [Filename; "file not found"]; %#ok<AGROW> 
                end
            end
            assert(length(Filename) == size(T0,1))
            T = addvars(T0, Filename, 'Before', 1);
            T(strcmp(T.Filename, "file not found"),:) = [];
        end
        function create_nifti_files_csv()
            cd(fullfile(getenv('ADNI_HOME'), 'bids', 'derivatives'))
            labels = {'cn', 'preclinical', 'cdr_0p5_aneg', 'cdr_0p5_apos', 'cdr_0p5_anan', 'cdr_gt_0p5_apos'};
            for idx = 1:length(labels)
                csv = sprintf('table_%s.csv', labels{idx});
                this = mladni.Anticlust(csv, label=labels{idx});
                this.writetable();
            end
        end
        function prepare_folders_for_VolBin()
            %% First run "anticlust_cn_repeat.R" created 20230526, modified 20230612

            pwd0 = pushd(fullfile(getenv("ADNI_HOME"), "bids", "derivatives"));
            g = glob("table_cn_rep*.csv");
            assert(length(g) == 100)
            for g1 = g'
                re = regexp(g1{1}, "table_cn_(?<rep>rep(A|B)\d+).csv", "names");
                fold = fullfile(getenv("ADNI_HOME"), "NMF_FDG", sprintf("baseline_cn_%s", re.rep));
                mkdir(fold);
                movefile(g1{1}, fullfile(fold, "nifti_files.csv"));
            end
            popd(pwd0);
        end
        function T = table_Filename_strrep(T, filename, str1, str2)
            %% TABLE_FILENAME_STRREP
            %  Args:
            %     T = []
            %     filename {mustBeTextScalar} = ''
            %     str1 {mustBeTextScalar} = '/home/usr/jjlee/mnt/CHPC_scratch'
            %     str2 {mustBeTextScalar} = '/scratch/jjlee'
            %  Returns:
            %     T

            arguments
                T = []
                filename {mustBeTextScalar} = ''
                str1 {mustBeTextScalar} = '/home/usr/jjlee/mnt/CHPC_scratch'
                str2 {mustBeTextScalar} = '/scratch/jjlee'
            end

            if isempty(T) && isfile(filename)
                T = readtable(filename, ReadVariableNames=false, Delimiter=',');
            end
            if isempty(T) && ~isfile(filename)
                return
            end

            assert(istable(T))
            if ~contains(T.Properties.VariableNames, 'Filename')
                T.Properties.VariableNames{1} = 'Filename';
            end
            Filename = strrep(T.Filename, str1, str2);
            if ~isempty(filename)
                writetable(table(Filename), filename, 'WriteVariableNames', false)
            end            
        end
    end

    methods
        function nums = num_each_category(~, var, cats)
            assert(iscell(var));
            assert(iscell(cats));
            N = numel(cats);
            for n = 1:N
                nums(n) = sum(strcmp(var, cats{n})); %#ok<AGROW> 
            end
        end
        function summary = summarise(this)
            %sumA = struct('MergeAge', struct('mean', mean(repA.MergeAge), 'std', std(repA.MergeAge)), 'MergePtGender', struct('Male', sum(strcmp(repA.MergePtGender, 'Male')), 'Female', sum(strcmp(repA.MergePtGender, 'Female'))), 'MergeDx', struct('CN', , 'MCI', 'Dementia', ))
            
            age = struct('repA', struct('mean', mean(this.repA.MergeAge), 'std', std(this.repA.MergeAge)), ...
                         'repB', struct('mean', mean(this.repB.MergeAge), 'std', std(this.repB.MergeAge)));

            num_A = this.num_each_category(this.repA.MergePtGender, {'Male', 'Female'});
            num_B = this.num_each_category(this.repB.MergePtGender, {'Male', 'Female'});
            sex = struct('repA', struct('num_M', num_A(1), 'num_F', num_A(2)), ...
                         'repB', struct('num_M', num_B(1), 'num_F', num_B(2)));

            num_A = this.num_each_category(this.repA.MergeDx, {'CN', 'MCI', 'Dementia'});
            num_B = this.num_each_category(this.repB.MergeDx, {'CN', 'MCI', 'Dementia'});
            dx = struct('repA', struct('num_CN', num_A(1), 'num_MCI', num_A(2), 'num_Dementia', num_A(3)), ...
                        'repB', struct('num_CN', num_B(1), 'num_MCI', num_B(2), 'num_Dementia', num_B(3)));

            summary = struct('MergeAge', age, 'MergePtGender', sex, 'MergeDx', dx);
        end        
        function writetable(this, varargin)
            ip = inputParser;
            addParameter(ip, "repA_fn",  fullfile(this.workpath, this.folder_repA, "nifti_files.csv"),  @istext);
            addParameter(ip, "repA_fn2", fullfile(this.workpath, this.folder_repA, "nifti_files2.csv"), @istext);
            addParameter(ip, "repB_fn",  fullfile(this.workpath, this.folder_repB, "nifti_files.csv"),  @istext);
            addParameter(ip, "repB_fn2", fullfile(this.workpath, this.folder_repB, "nifti_files2.csv"), @istext);
            parse(ip, varargin{:});
            ipr = ip.Results;

            ensuredir(fileparts(ipr.repA_fn));
            ensuredir(fileparts(ipr.repB_fn));

            repA_flist  = strrep(this.repA.Filename, "/mnt/CHPC_scratch", "/scratch/jjlee");
            repA_flist2 = strrep(this.repA.Filename, "/mnt/CHPC_scratch", "/home/usr/jjlee/mnt/CHPC_scratch");
            repB_flist  = strrep(this.repB.Filename, "/mnt/CHPC_scratch", "/scratch/jjlee");
            repB_flist2 = strrep(this.repA.Filename, "/mnt/CHPC_scratch", "/home/usr/jjlee/mnt/CHPC_scratch");

            writetable(table(repA_flist),  ipr.repA_fn,  'WriteVariableNames', false);
            writetable(table(repA_flist2), ipr.repA_fn2, 'WriteVariableNames', false);
            writetable(table(repB_flist),  ipr.repB_fn,  'WriteVariableNames', false);
            writetable(table(repB_flist2), ipr.repB_fn2, 'WriteVariableNames', false);
        end

        function this = Anticlust(varargin)
            %% ANTICLUST 
            %  Args:
            %      original_csv (optional file): csv containing mladni.FDGDemographics.table_covariates; 
            %                                    default := "anticlust.csv"
            
            ip = inputParser;
            addOptional(ip, "original_csv", "anticlust.csv", @isfile);
            addParameter(ip, "repA_csv", "", @istext);
            addParameter(ip, "repB_csv", "", @istext);
            addParameter(ip, "label", "", @istext);
            addParameter(ip, "workpath", fullfile(getenv('ADNI_HOME'), 'NMF_FDG'), @isfolder)
            parse(ip, varargin{:})
            ipr = ip.Results;
            if ipr.repA_csv == ""
                ipr.repA_csv = strcat(myfileprefix(ipr.original_csv), "_repA.csv");
            end
            if ipr.repB_csv == ""
                ipr.repB_csv = strcat(myfileprefix(ipr.original_csv), "_repB.csv");
            end
            this.label = ipr.label;
            this.workpath = ipr.workpath;
            
            this.repA = readtable(ipr.repA_csv, 'VariableNamesLine', 1);
            this.repB = readtable(ipr.repB_csv, 'VariableNamesLine', 1);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
