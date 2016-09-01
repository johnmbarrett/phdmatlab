function mungeWaveClusToNex(recording,channels,suffix)
    if nargin < 2
        channels = [];
    end
    
    [~,filename] = getAnalysisOutputDir(recording);
    
    nexFile = nexCreateFileData(25000);
    
    function fn(spikeTimes,~,channelLabel,~,~,spikes)
        nexFile = nexAddWaveform(nexFile,25000,spikeTimes,1e3*spikes',['ch' channelLabel]);
    end
    
    forEachChannel(recording,channels,false,@fn,-1,true);
    
    if nargin > 2
        fileSuffix = suffix;
    elseif isempty(channels)
        fileSuffix = '';
    elseif isnumeric(channels) && isequal(channels,channels(1):channels(end))
        fileSuffix = ['_chs_' num2str(channels(1)) '-' num2str(channels(end))];
    else
        fileSuffix = 'chs';
        
        if ischar(channels)
            n = size(channels,1);
            indexFn = @(x,ii) x(ii,:);
        else
            n = numel(channels);
            
            if isnumeric(channels)
                indexFn = @(x,ii) num2str(x(ii));
            elseif iscell(channels)
                indexFn = @(x,ii) x{ii};
            else
                error('Unsupported channel specification');
            end
        end
        
        for ii = 1:n
            fileSuffix = [fileSuffix '_' indexFn(channels,ii)]; %#ok<AGROW>
        end
    end
    
    writeNexFile(nexFile,[filename ' spikes trimmed ' fileSuffix '.nex']);
end