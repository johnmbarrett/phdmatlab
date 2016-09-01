function [data,sampleCount,eof] = getMCSRawData(file,startIndex,indexCount,channels,timeUnits,origin)
    if nargin < 3
        error('You must provide a file structure, start time and length of data segment');
    end
    
    if nargin < 6
        origin = 0;
    elseif ~isnumeric(origin) || ~isscalar(origin) || ~ismember(origin,[-1 0 1])
        error('Origin argument must be -1 (beginning of file), 0 (current position), or 1 (end of file)');
    end
    
    if nargin >= 5 % default to samples
        if ischar(timeUnits)
            if strncmpi(timeUnits,'se',2)
                % startIndex in [0,1/sampleRate) corresponds to the first
                % sample, which is index 0 (unlike most Matlab arrays)
                startIndex = floor(startIndex*file.sampleRate); 
                % indexCount in [0,1/sampleRate) is asking for 0 samples
                indexCount = floor(indexCount*file.sampleRate);
            elseif strncmpi(timeUnits,'sa',2)
                startIndex = startIndex-1; % Matlab arrays are indexed from 1, files from 0
            else
                error('Invalid unit specifier %s\n',timeUnits);
            end
        else
            error('Units argument must be a string beginning ''se'' (for seconds) or ''sa'' (for samples)');
        end
    end
    
    assert(indexCount > 0,'Data segment must be longer than zero samples (%d samples requested)\n',indexCount);
    
    labels = cellstr(vertcat(file.electrodes.label));
    
    % TODO : other types of stream?
    if nargin < 4 || isempty(channels)
        channels = vertcat(file.electrodes.index);
    elseif isnumeric(channels)
        assert(all(isfinite(channels)) && all(channels > 0) && all(channels <= file.nChannels),'Numeric channel specifiers must be within the range of valid channels indices (1 to %d)\n',file.nChannels);
    elseif ischar(channels)
        channelSpec = channels;
        channels = [];
        
        for ii = 1:size(channelSpec,1)
            channelIndex = find(strcmpi(channelSpec(ii,:),labels));
            
            if isempty(channelIndex)
                warning('Unable to locate channel %s\n',channelSpec(ii,:));
                continue;
            elseif numel(channelIndex) > 1
                error('Multiple channel matches for channel %s, file may be corrupt\n',channelSpec(ii,:));
            end
            
            channels(end+1) = file.electrodes(channelIndex).index; %#ok<AGROW>
        end
    elseif iscell(channels)
        channelSpec = channels;
        channels = [];
        
        for ii = 1:size(channelSpec,1)
            % TODO : code duplication; refactor
            if isnumeric(channelSpec{ii})
                channelIndex = channelSpec{ii};
                assert(all(isfinite(channelIndex)) && all(channelIndex > 0) && all(channelIndex <= file.nChannels),'Numeric channel specifiers must be within the range of valid channels indices (1 to %d)\n',file.nChannels);
                channels = [channels channelIndex]; %#ok<AGROW>
            elseif ischar(channelSpec{ii})
                channelIndex = find(strcmpi(channelSpec{ii},labels));
                
                if isempty(channelIndex)
                    warning('Unable to locate channel %s\n',channelSpec{ii});
                    continue;
                elseif numel(channelIndex) > 1
                    error('Multiple channel matches for channel %s, file may be corrupt\n',channelSpec{ii});
                end

                channels(end+1) = file.electrodes(channelIndex).index; %#ok<AGROW>
            end
        end
    end
    
    channels = channels(:); % probably not important but feels neat
    
    status = fseek(file.handle,startIndex*file.nStreams*2+file.headerSize*(origin == -1),origin);
    
    if status < 0
        error('Unexpected error seeking to sample #%d\n',startIndex+1);
    elseif feof(file.handle)
        error('Unexpected EOF seeking to sample #%d\n',startIndex+1);
    end
    
    [data,count] = fread(file.handle,[file.nStreams indexCount],'int16');
    eof = feof(file.handle);
    
    sampleCount = ceil(count/file.nStreams);
    
    if sampleCount ~= count/file.nStreams
        warning('Unable to fetch the requested range for all channels, file may be corrupt');
    end
    
    data = data(channels,:)'*file.scaling;
end