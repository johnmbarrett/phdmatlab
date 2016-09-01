function oliverSpikesToNex2(hdf5file,outputFile,maxChannels)
    if ~ismember(nargin,[2 3])
        fprintf('Usage: oliverSpikesToNex2([extracted spikes file],[output file name],sampling frequency in Hz,max # channels\n');
        return;
    end
    
    if nargin == 2 || numel(maxChannels) ~= 1 || ~isnumeric(maxChannels) || isnan(maxChannels) || maxChannels < 1
        maxChannels = Inf;
    end

    if isempty(outputFile)
        outputFile = uiputfile('*.nex','Choose .nex file to save to...');
    end
    
    outputFile = regexprep(outputFile,'.nex$','');

    wb = waitbar(0,'Loading spikes...');
    
    if isempty(hdf5file)
        hdf5file = uigetfile('*.hdf5','Choose file with extracted spikes...');
    elseif ~exist(hdf5file,'file')
        error('Could not find spikes file: %s',shapeFile);
    end
    
    waitbar(0,wb,'Loading extracted spikes...');
        
    try
        channels = h5read(hdf5file,'/Channels');
        waitbar(1/4,wb);
        times = h5read(hdf5file,'/Times');
        waitbar(2/4,wb);
        shapes = h5read(hdf5file,'/Shapes');
        waitbar(3/4,wb);
        sampleRate = h5read('P26_23May13_DA_WN.hdf5','/Sampling');
        waitbar(4,wb);
    catch err
        logMatlabError(err);
        error 'Could not extract spikes from selected hdf5 file; please choose another one.';
    end
    
    uniqueChannels = unique(channels);
    nChannels = numel(uniqueChannels);
    
    nexFile = nexCreateFileData(sampleRate);
    
    nFiles = 1;
    processedChannels = 0;
    
    for ii = 1:nChannels
        waitbar((ii-1)/nChannels,wb,['Processing channel #' num2str(ii) ' of ' num2str(nChannels)]);
        
        channel = uniqueChannels(ii);
        channelSpikes = channels == channel;
        channelTimes = times(channelSpikes)/sampleRate;
        channelShapes = shapes(:,channelSpikes);
        
        channelLabel = sprintf('Ch%02d_%02d',mod(channel,64)+1,floor(channel/64)+1);
        nexFile = nexAddWaveform(nexFile,sampleRate,channelTimes(:),channelShapes,channelLabel);
        
        processedChannels = processedChannels + 1;
            
        if processedChannels == maxChannels || ii == nChannels
            writeNexFile(nexFile,sprintf('%s_%02d.nex',outputFile,nFiles));
            processedChannels = 0;
            nFiles = nFiles + 1;
            nexFile = nexCreateFileData(sampleRate);
        end
    end
    
    close(wb);
end