classdef Anticlust
    %% line1
    %  line2
    %  
    %  Created 18-Aug-2022 21:27:03 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2009381 (R2022a) Update 4 for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        repA
        repB
    end

    methods
        function call(this)
        end
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
            addParameter(ip, "repA_fn",  "repA_nifti_files.csv",  @istext);
            addParameter(ip, "repA_fn2", "repA_nifti_files2.csv", @istext);
            addParameter(ip, "repB_fn",  "repB_nifti_files.csv",  @istext);
            addParameter(ip, "repB_fn2", "repB_nifti_files2.csv", @istext);
            parse(ip, varargin{:});
            ipr = ip.Results;

            repA_flist  = strrep(this.repA.Filelist, "/mnt/CHPC_scratch", "/scratch/jjlee");
            repA_flist2 = this.repA.Filelist;
            repB_flist  = strrep(this.repB.Filelist, "/mnt/CHPC_scratch", "/scratch/jjlee");
            repB_flist2 = this.repB.Filelist;

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
            parse(ip, varargin{:})
            ipr = ip.Results;
            if ipr.repA_csv == ""
                ipr.repA_csv = strcat(myfileprefix(ipr.original_csv), "_repA.csv");
            end
            if ipr.repB_csv == ""
                ipr.repB_csv = strcat(myfileprefix(ipr.original_csv), "_repB.csv");
            end
            
            this.repA = readtable(ipr.repA_csv);
            this.repB = readtable(ipr.repB_csv);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
