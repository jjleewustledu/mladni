classdef FDGQC < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% line1
    %  line2
    %  
    %  Created 21-Mar-2022 14:24:57 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.10.0.1851785 (R2021a) Update 6 for MACI64.  Copyright 2022 John J. Lee.
    methods (Static)
        function propcluster()            
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '4000';
            c.AdditionalProperties.Node = 10;
            c.AdditionalProperties.Partition = 'test';
            c.AdditionalProperties.WallTime = '1:00:00';
            c.saveProfile
        end
        function getDebugLog(j,c)
            try
                c.getDebugLog(j)
            catch
                c.getDebugLog(j.Tasks(end))
            end
        end
        function [j,c] = parcluster()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            
            c = parcluster;
            disp(c.AdditionalProperties)
            j = c.batch(@mladni.FDGQC.batch, 1, {}, 'Pool', 20, ...
                'CurrentFolder', '.', 'AutoAddClientPath', false);
            
        end
        function t = batch(varargin)
            %% #BATCH
            
            ip = inputParser;
            addParameter(ip, 'proc', 'CASU', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mladni.FDG.batch()')        
            
            t0 = tic;            
            setenv('ADNI_HOME', '/home/aris_data/ADNI_FDG')
            globbed = globT( ...
                fullfile( ...
                    '/home/aris_data/ADNI_FDG/bids/derivatives', ...
                    'sub-*', 'ses-*', 'pet', ...
                    sprintf('*-%s-orientstd_pet_final.json', ipr.proc)));
            fprintf('mladni.FDG.batch.globbed.size:\n')
            disp(size(globbed))
            rerr = nan(size(globbed));
            terr = nan(size(globbed));
            fcost = nan(size(globbed));
            
            parfor idx = 1:length(globbed)
                try
                    j = jsondecode(fileread(globbed{idx}));
                    rerr(idx) = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_rotation_error.err); %#ok<PFOUS>
                    terr(idx) = str2double(j.mlfourdfp_SimpleT4ResolveBuilder.cost_final.pairs_translation_error.err); %#ok<PFOUS>
                    fcost(idx) = j.mlfsl_Flirt.cost_final.cost; %#ok<PFOUS>
                catch ME
                    handwarning(ME);
                end
            end
            save(fullfile(getenv('ADNI_HOME'), 'rerr.mat'), 'rerr')
            save(fullfile(getenv('ADNI_HOME'), 'terr.mat'), 'terr')
            save(fullfile(getenv('ADNI_HOME'), 'fcost.mat'), 'fcost')
            t = toc(t0);

            disp('mladni.FDG.batch() completed')            
        end
    end
    
    methods
        function this = FDGQC(varargin)
            %% FDGQC 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            ip = inputParser;
            addParameter(ip, "arg1", [], @(x) false)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
