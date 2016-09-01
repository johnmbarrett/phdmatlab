function oliverSpikesToNex(matlabSpikeFile,asciiSpikeFile,shapeFile,outputFile,sampleRate,maxChannels)
    if nargin ~= 6
        fprintf('Usage: oliverSpikesToNex([extracted spikes file],[original spikes file],[spike shapes file],[output file name],sampling frequency in Hz\n');
        return;
    end

    if isempty(outputFile)
        outputFile = uigetfile('*.nex','Choose .nex file to save to...');
    end
    
    outputFile = regexprep(outputFile,'.nex$','');

    wb = waitbar(0,'Loading spike shapes...');
    
    if isempty(shapeFile)
        shapeFile = uigetfile('*.txt','Choose file with spike shapes...');
    elseif ~exist(shapeFile,'file')
        error('Could not find shape file: %s',shapeFile);
    end
    
    shapes = dlmread(shapeFile);
    
    waitbar(1/3,wb,'Loading extracted spike times...');
    
    if isempty(matlabSpikeFile)
        matlabSpikeFile = uigetfile('*.mat','Choose file with extracted spike times...');
    elseif ~exist(matlabSpikeFile,'file')
        error('Could not find Matlab spike file: %s',matlabSpikeFile);
    end
    
    load(matlabSpikeFile);
    
    waitbar(2/3,wb,'Loading original spike times...');
    
    if isempty(asciiSpikeFile)
        asciiSpikeFile = uigetfile('*.txt','Choose file with original spike times...');
    elseif ~exist(asciiSpikeFile,'file')
        error('Could not find ASCII spike file: %s',asciiSpikeFile);
    end
    
    spikes = dlmread(asciiSpikeFile);
    
    close(wb);
    
    nexFile = nexCreateFileData(sampleRate);
    
    msgFn = @(x) ['Processing channel #' num2str(x) ' of 4096'];
    wb = waitbar(0,msgFn(0));
    nChannels = 0;
    nFiles = 1;
    for ii = 1:64
        for jj = 1:64
            channelIndex = (ii-1)*64+jj-1;
            waitbar((channelIndex+1)/4096,wb,msgFn(channelIndex + 1));
            
            spikeVar = ['Ch' num2str(ii) '_' num2str(jj)];
            
            if ~exist(spikeVar,'var')
                continue;
            end
            
            eval(['spikeTrain = ' spikeVar ';']);
            spikeTimes = find(spikeTrain ~= 0);
            
            if isempty(spikeTimes)
                continue;
            end
            
            nChannels = nChannels + 1;
            
            spikeIndices = (spikes(:,1) == channelIndex) & (ismember(spikes(:,2),spikeTimes));
            waveforms = shapes(spikeIndices,:);
            
            if size(waveforms,1) ~= numel(spikeTimes)
                if size(waveforms,1) < numel(spikeTimes)
                    problem = 'Missing';
                else
                    problem = 'Extra';
                end
                
                warning('%s spikes on channel (%d,%d) in file %s: expected %d spikes, found %d spikes',problem,ii,jj,matlabSpikeFile,numel(spikeTimes),size(waveforms,1)); %#ok<WNTAG>
            end
            
            nexFile = nexAddWaveForm(nexFile,sampleRate,full(spikeTimes/sampleRate),waveforms'/1e3,spikeVar);
            
            if nChannels > maxChannels
                writeNexFile(nexFile,sprintf('%s_%02d.nex',outputFile,nFiles));
                nFiles = nFiles + 1;
                nChannels = 0;
                nexFile = nexCreateFileData(sampleRate);
            end
        end
    end
    close(wb);
    writeNexFile(nexFile,sprintf('%s_%02d.nex',outputFile,nFiles));
end