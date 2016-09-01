function [indicess,channels,infos] = getMCDChannels(file,fileInfo,type,varargin)
    opt = getopt('channels=NaN chlabels=NaN',varargin{:});

    if iscell(type)
        validIndices = [];
        
        for ii = 1:numel(type)
            if ~ischar(type{ii})
                warning('Type specifier %d is not a valid identity specifier and will be ignored:\n',ii);
                disp(type{ii});
            else
                validIndices(end+1) = ii; %#ok<AGROW>
            end
        end
        
        type = type(validIndices);
        indicess = cell(size(type));
        channels = cell(size(type));
        infos = cell(size(type));
    elseif ischar(type)
        type = {type};
        indicess = {[]};
        channels = {[]};
        infos = {[]};
    end

    [err,entities] = ns_GetEntityInfo(file,1:fileInfo.EntityCount);
    
    if err
        error('Unable to get entity info'); % TODO : name???
    end
    
    for ii = 1:numel(type)
        specifier = type{ii}(1);
        
        filter = @(x) true(size(x));
        if strncmpi(specifier,'a',1)
            typeNumber = 2;
            fn = @ns_GetAnalogInfo;
        elseif strncmpi(specifier,'el',1)
            typeNumber = 2;
            fn = @ns_GetAnalogInfo;
            filter = @(entities) strncmpi({entities.EntityLabel}','elec',4);
        elseif strncmpi(specifier,'ev',1)
            typeNumber = 1;
            fn = @ns_GetEventInfo;
        elseif strncmpi(specifier,'n',1)
            typeNumber = 4;
            fn = @ns_GetNeuralInfo;
        elseif strncmpi(specifier,'s',1)
            typeNumber = 3;
            fn = @ns_GetSegmentInfo;
        elseif strncmpi(specifier,'u',1)
            typeNumber = 0;
            fn = @(file,entityIDs) struct(zeros(numel(entityIDs,0)));
        else
            warning('Invalid type specifier %s\n',specifier);
            channels{ii} = [];
            continue;
        end
        
        indices = find(vertcat(entities.EntityType) == typeNumber & filter(entities));
        
        if isnumeric(opt.channels) && ~isnan(opt.channels) && isvector(opt.channels) && min(opt.channels) >= 1 && max(opt.channels) <= numel(indices)
            indices = indices(opt.channels);
        end
        
        if iscell(opt.chlabels)
            isSelected = false(size(indices));
            
            for jj = 1:numel(indices)
                isSelected(jj) = any(strcmp(entities(indices(jj)).EntityLabel(end-1:end),opt.chlabels));
            end
               
            indices = indices(isSelected);
        end
            
        indicess{ii} = indices;
        channels{ii} = entities(indices);
        [~,infos{ii}] = fn(file,indices);
    end
end