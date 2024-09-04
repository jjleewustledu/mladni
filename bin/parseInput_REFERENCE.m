function param = parseInput(mfMeth,isList,numBases,varargin)

% function returns a param structure with the following fields
%
%   mfMeth   : determines which matrix factorization method is going to be
%              used
%   isList   : determines if the input is given in a .lst file or as a matlab
%   numBases : determines the number of bases
%   w0       : initialization of W
%   initMeth : methods that specifies how w0 should be initializad
%   seed     : seeds random initialization
%   ouputdir : directory where to write the results
%   forseSave: 0 or 1, forses the saving of the results
%   saveInterm: 0 or 1, saves intermediate results
%   sparse   : [0, 1) where 0 means completely distributed and 1 means 
%              ultimate sparseness.
% downSample: factor by which data are downsampled (default 1 : no
%             downsampe)
%   permute : determines if the variables of each sample are going to be
%             permuted or not (0 or 1)
%   permSeed : seed used for random permutation
%   idx      : indices that determine which sample should considered during
%              the training (used when performing Cross Validation). The
%              indices are given in a .txt file that has been generated in
%              the output directory of the current fold experiment
%   dMapDist : It specifies how to calculate the distance for the diffusion 
%              maps. This is used in the case that KMeans with Diffusion
%              Maps are used
%   dMapSup  : 0 or 1, it specifies if label information is taken into 
%              account when calculating the distances
%   negPos  : specifies how to handle negative values in the data. It can be
%             either 0, 1, or 2. If 0, then negative values are thresholded out
%             (negative values = 0). If 1, then negative values become
%             positive by adding appropriate value (default). If 2, then
%             negative values become positive by taking absolute value
%   smooth  : specifies the kernel that is going to be used to smooth the
%             input data

% initialize parser
p = inputParser;
p.addParamValue('w0',[],@(x) true);
p.addParamValue('initMeth',0,@(x) true);
p.addParamValue('seed',[],@(x) true);
p.addParamValue('outputDir','./',@(x) true);
p.addParamValue('forseSave',0,@(x) true);
p.addParamValue('saveInterm',0,@(x) true);
p.addParamValue('sparse',[],@(x) true);
p.addParamValue('downSample',1,@(x) true);
p.addParamValue('permute',0,@(x) true);
p.addParamValue('permSeed',1,@(x) true);
p.addParamValue('idx',[],@(x) true);
p.addParamValue('dMapDist',[],@(x) true);
p.addParamValue('dMapSup',0,@(x) true);
p.addParamValue('negPos',1,@(x) true);
p.addParamValue('smooth',0,@(x) true);
p.addParamValue('mask',[],@(x) true);

dflts = {[] 0 [] './' 0 0 [] 1 0 1 [] [] 0 1 0 []} ;


% parse input depending on whether code is running in deployed or MATLAB mode 
if(isdeployed)        
    % obligatory input arguments
    param.mfMeth = mfMeth;
    param.isList = str2double(isList) ;
    disp(['isList ' num2str(param.isList)]);    
    param.numBases = str2double(numBases);  %numBases is integer
    disp(['numBases ' num2str(param.numBases)]);
    
    % parse input parameters
    if (~isempty(varargin))
        p.parse(varargin{:});        
        % if optional parameters are given       
        if(any(strcmp('initMeth',p.UsingDefaults)))
            param.initMeth = p.Results.initMeth ;
            disp(['Default initMeth value: ' num2str(param.initMeth)]);
        else
            param.initMeth = str2double(p.Results.initMeth);
            disp(['User initMeth value: ' num2str(param.initMeth)]);
        end        
        if(any(strcmp('seed',p.UsingDefaults)))            
            param.seed = p.Results.seed ;
            disp('Default seed value: empty (shuffle option for random generator)');
        else
            param.seed = str2double(p.Results.seed);
            disp(['User seed value: ' num2str(param.seed)]);            
        end
        if(any(strcmp('permSeed',p.UsingDefaults)))            
            param.permSeed = p.Results.permSeed ;
            disp('Default seed value for permutation: 0');
        else
            param.permSeed = str2double(p.Results.permSeed);
            disp(['User seed value for permutation: ' num2str(param.permSeed)]);            
        end
        if(any(strcmp('sparse',p.UsingDefaults)))
            param.sparse = p.Results.sparse ;
            disp('Default sparse value: empty');
        else
            param.sparse = str2double(p.Results.sparse);
            disp(['User seed value: ' num2str(param.sparse)]);
        end
        if(any(strcmp('forseSave',p.UsingDefaults)))            
            param.forseSave = 1 ;
            disp('Deployed case: results will be saved!');   
        else
            param.forseSave = 1 ;
            disp('Deployed case: results will be saved!')
        end
        if(any(strcmp('saveInterm',p.UsingDefaults)))            
            param.saveInterm = p.Results.saveInterm ;
            disp('Default saveInterm value : no intermediate results will be saved');          
        else
            param.saveInterm = str2double(p.Results.saveInterm) ;
            disp(['User saveInterm flag: ' num2str(param.saveInterm)]) ;
        end        
        if(any(strcmp('negPos',p.UsingDefaults)))
            param.negPos = p.Results.negPos ;
            disp('Default negPos value : negative values will become positive by adding appropriate value');
        else
            param.negPos = str2double(p.Results.negPos) ;
            disp(['User negPos flag: ' num2str(param.negPos)]) ;
        end
        if(any(strcmp('permute',p.UsingDefaults)))            
            param.permute = p.Results.permute ;
            disp('Default permute value : no data permutation will be performed');          
        else
            param.permute = str2double(p.Results.permute) ;
            disp(['User permute flag: ' num2str(param.permute)]) ;
        end     
        if(any(strcmp('downSample',p.UsingDefaults)))
            param.downSample = p.Results.downSample ;
            disp(['Default downSample value: ' num2str(param.downSample)]);
        else
            param.downSample = str2double(p.Results.downSample);
            disp(['User downSample value: ' num2str(param.downSample)]);
        end
        if(any(strcmp('dMapSup',p.UsingDefaults)))
            param.dMapSup = p.Results.dMapSup ;
            disp('Default dMapSup value: 0 - no supervision');
        else
            param.dMapSup = str2double(p.Results.dMapSup);
            disp(['User dMapSup value: ' num2str(param.dMapSup)]);
        end
        if(any(strcmp('smooth',p.UsingDefaults)))
            param.smooth = p.Results.smooth ;
            disp('Default smooth value: 0 - no smoothing will be performed');
        else
            param.smooth = str2double(p.Results.smooth);
            disp(['User smooth value: ' num2str(param.smooth)]);
        end
        param.idx = p.Results.idx ;
        param.outputdir = p.Results.outputDir ;
        param.mask = p.Results.mask ;
        param.dMapDist  = p.Results.dMapDist ;
        if(~isempty(p.Results.w0))
            tmp_struct = load(p.Results.w0,'w0');
            param.w0 = tmp_struct.w0 ; clear tmp_struct
        else
            param.w0 = p.Results.w0 ; % empty
        end
        clear p ; 
        % varargin is empty and default values are used
    else
        [param.w0, param.initMeth, param.seed, param.outputdir param.forseSave, param.saveInterm ...
            param.sparse param.downSample param.permute param.permSeed param.idx param.dMapDist ...
            param.dMapSup param.negPos param.smooth param.mask] = dflts{:} ;
        param.forseSave = 1 ;
    end % end if(nargin>4)
    % when deployed, an output directory should be given
    if(~isempty(param.outputdir))
        disp(['Output directory: ' param.outputdir]);
    else
        error('extractBasesPNMF:argChk','a directory to save the results should be provided!');
    end
else
    param.mfMeth = mfMeth ;
    param.isList = isList ;
    param.numBases = numBases ;
    disp(['isList ' num2str(param.isList)]);
    disp(['numBases ' num2str(param.numBases)]);
    % In the case that is used as simple matlab function, the folder should
    % be added to the path, or:
    if (~isempty(varargin))
        p.parse(varargin{:});
        param.initMeth   = p.Results.initMeth ;
        param.seed       = p.Results.seed ;
        param.w0         = p.Results.w0 ;
        param.outputdir  = p.Results.outputDir ;
        param.forseSave  = p.Results.forseSave ;
        param.saveInterm = p.Results.saveInterm ;
        param.sparse     = p.Results.sparse ;
        param.downSample = p.Results.downSample ;
        param.permute    = p.Results.permute ;
        param.permSeed   = p.Results.permSeed ;
        param.idx        = p.Results.idx ;
        param.dMapDist   = p.Results.dMapDist ;
        param.dMapSup    = p.Results.dMapSup ;
        param.negPos     = p.Results.negPos ;
        param.smooth     = p.Results.smooth ;
        param.mask       = p.Results.mask ;
        % varargin is empty and default values are used
    else
        [param.w0, param.initMeth, param.seed, param.outputdir, param.forseSave, param.saveInterm, ...
            param.sparse, param.downSample, param.idx, param.permute, param.permSeed, param.dMapDist, ...
            param.dMapSup, param.negPos, param.smooth, param.mask] = dflts{:} ;
    end
    clear p;
end % end if(isdeployed)
