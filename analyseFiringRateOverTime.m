function analyseFiringRateOverTime(recordings,channels,varargin)
    nRecordings = numel(recordings);
    nChannels = numel(channels);
    
    nSpikes = zeros(nChannels,nRecordings);
    
    function fn(spikeTimes,channelIndex,varargin)
        nSpikes(channelIndex,ii) = numel(spikeTimes);
    end
    
    for ii = 1:nRecordings
        recording = recordings{ii};
        
        forEachChannel(recording,channels,true,@fn,false)
    end
    
    function missingAllSpikesFilesFn(spikeTimes,channelIndex,~,cluster,varargin)
        if cluster == 0
            return;
        end
        
        if size(nSpikes,1) < channelIndex
            nSpikes(channelIndex,ii) = numel(spikeTimes);
        else
            nSpikes(channelIndex,ii) = nSpikes(channelIndex,ii) + numel(spikeTimes);
        end
    end

    if isempty(nSpikes) || sum(sum(abs(nSpikes))) == 0
        for ii = 1:nRecordings
            recording = recordings{ii};

            forEachChannel(recording,channels,true,@missingAllSpikesFilesFn,true)
        end
    end
    
    currentDir = pwd;
    slashes = strfind(currentDir,'\');
    
    if ~isempty(slashes)
        currentDir = currentDir(slashes(end)+1:end);
    end
    
    figure;
    boxplot(nSpikes,recordings(:));
    xlabel('Recording');
    ylabel('# Spikes');
    title(currentDir);
    
    options = getopt('suffix=NaN',varargin{:});
    
    suffix = '';
    if ischar(options.suffix)
        suffix = ['_' options.suffix];
    end
    
    figFile = sprintf('%s%s_firing_rates',currentDir,suffix);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    
    spreadsheet = [figFile '.xlsx'];
    xlswrite(spreadsheet,recordings,sprintf('A1:%s1','A'+numel(recordings)-1));
    xlswrite(spreadsheet,nSpikes,sprintf('A2:%s%d','A'+numel(recordings)-1,size(nSpikes,1)+1));
end