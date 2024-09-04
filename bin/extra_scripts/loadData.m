function data = loadData(dataInput,param,subsetIdx)
     
if(param.isList==1)
    disp(['List ' dataInput ]);    
    data = loadDataFromList(dataInput,param,subsetIdx) ; % data is a struct !
else
    % if data input is not given as list, then what we do depends on
    % whether the code is used as a function or a standalone executable
    if(isdeployed)
        % dataInput contains the full path to the .mat that contains
        % data
        tmp = load(dataInput,'data');
        tmp = struct2cell(tmp);
        data = tmp{1}; clear tmp
        data.dsflag = 0;                
    else
        % dataInput is the data_matrix
        data = dataInput ; clear dataInput
    end
    % we have loaded data but we don't know if it is a struct or just an
    % array
    if (isstruct(data))
        % check if field X is there
        if(~isfield(data,'X'))
            error('extractBases:argChk','The field X that should contain the data is missing !');
        end
        if(~isempty(subsetIdx))
            data.X = data.X(:,subsetIdx==1);
        end
        if(isfield(data,'dimx') && isfield(data,'dimy') && isfield(data,'dimz') && (param.downSample~=1))
            % downsample data
            tmp = data.X ;
            h_nsamples = size(data.X,2) ;
            data.X = [] ;
            new_size_x = ceil(data.dimx/param.downSample) ;
            new_size_y = ceil(data.dimy/param.downSample) ;
            new_size_z = ceil(data.dimz/param.downSample) ;
            data.X = zeros(new_size_x*new_size_y*new_size_z,h_nsamples) ;
            for i=1:h_nsamples
                xx1 = linspace(1,data.dimx,data.dimx) ; 
                yy1 = linspace(1,data.dimy,data.dimy) ;
                zz1 = linspace(1,data.dimz,data.dimz) ;
                xx2 = linspace(1,data.dimx,new_size_x) ;
                yy2 = linspace(1,data.dimy,new_size_y) ;
                zz2 = linspace(1,data.dimz,new_size_z) ;
                [XX1,YY1,ZZ1] = meshgrid(yy1,xx1,zz1) ;
                [XX2,YY2,ZZ2] = meshgrid(yy2,xx2,zz2) ;
                h_img = reshape(tmp(:,i),data.dimx.data.dimy,data.dimz) ;
                img_interp = interp3(XX1,YY1,ZZ1,h_img,XX2,YY2,ZZ2) ;
                data.X(:,i) = double(img_interp(:));
            end
            clear tmp
            clear h_nsamples
            clear h_img
            data.dsflag = 1;
        else
            data.dsflag = 0;
        end
    else
        % put data into a struct
        tmp = data ; clear data ; data.X = tmp ; clear tmp ;
        if(~isempty(subsetIdx))
            data.X = data.X(:,subsetIdx==1);
        end
        data.dsflag = 0;
    end
end

if(param.permute)
    s = RandStream.create('mt19937ar','seed',param.permSeed);RandStream.setGlobalStream(s);
    for i=1:size(data.X,2)
        tmp = data.X(randperm(length(data.X(:,i))),i);
        data.X(:,i) = tmp ;
    end
end
