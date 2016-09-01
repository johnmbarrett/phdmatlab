function [out,ch,cl] = forEachChannel(recording,channelInfo,ignoreNoise,fn,useClustered,isTrimmed)
    if nargin < 6
        isTrimmed = true;
    end
    
    % useClustered is trinary: 1 = force clustered, 0 = clustered if
    % available, -1 = force unclustered.
    % TODO : the else statement here is backwards compatability for the
    % legacy syntax, assuming that if it's unspecified then we want
    % clustered if available, if it's specified false then we definitely
    % want unclustered and if it's specified true then we definitely want
    % clustered.  Should really go through all usages and check, but MATLAB
    % is not a good IDE.
    if nargin < 5
        useClustered = 0;
    elseif islogical(useClustered)
        useClustered = 2*useClustered-1;
    end
    
    if isempty(channelInfo)
        channelInfo = struct('label',cellstr(num2str(channelIndexToMCSChannelNumber(1:60)')));
    elseif isnumeric(channelInfo)
        channelInfo = struct('label',cellstr(num2str(channelIndexToMCSChannelNumber(channelInfo(:)))));
    elseif ischar(channelInfo)
        channelInfo = struct('label',cellstr(channelInfo));
    elseif iscell(channelInfo)
        channelInfo = struct('label',channelInfo);
    end
    
    if nargout > 3
        error no
    end
    
    if nargout == 3
        cl = [];
    end
    
    if nargout >= 2
        ch = {};
    end
    
    if nargout >= 1
        out = {};
    end
    
    for channelIndex = 1:numel(channelInfo)
        channelLabel = channelInfo(channelIndex).label;
        [spikeFile,isClustered] = getSpikeFile(recording,channelLabel,useClustered,isTrimmed);
        
        if ~exist(spikeFile,'file')
            continue;
        end
        
        if isClustered
            load(spikeFile,'cluster_class','spikes');
        else
            load(spikeFile,'index','spikes');
            cluster_class = [ones(numel(index),1) index(:)]; %#ok<NODEF>
        end
        
        clusters = unique(cluster_class(:,1));
        
        for ii = 1:numel(clusters)
            cluster = clusters(ii);
            
            if ignoreNoise && cluster == 0
                continue;
            end
            
            spikeTimes = cluster_class(cluster_class(:,1) == cluster,2)/100;
            
            if nargout == 3
                cl(end+1) = cluster; %#ok<AGROW>
            end
            
            if nargout >= 2
                ch{end+1} = channelLabel; %#ok<AGROW>
            end
            
            if nargout >= 1
                out{end+1} = fn(spikeTimes,channelIndex,channelLabel,cluster,max(clusters),spikes); %#ok<AGROW>
            else
                fn(spikeTimes,channelIndex,channelLabel,cluster,max(clusters),spikes);
            end
        end
    end
end