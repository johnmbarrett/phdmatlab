function matthiasSpikesToNex(matlabSpikeFile,asciiSpikeFile,shapeFile,outputFile,sampleRate)
    if ~exist(shapeFile,'file')
        error('Could not find shape file: %s',shapeFile);
    end
    
    shapes = dlmread(shapeFile);
    
    if ~exist(matlabSpikeFile,'file')
        error('Could not find Matlab spike file: %s',matlabSpikeFile);
    end
    
    load(matlabSpikeFile);
    
    if ~exist(asciiSpikeFile,'file')
        error('Could not find ASCII spike file: %s',asciiSpikeFile);
    end
    
    spikes = dlmread(asciiSpikeFile);
    
    nexFile = nexCreateFileData(sampleRate);
    
    for ii = 1:64
        for jj = 1:64
            spikeVar = ['Ch' num2str(ii) '_' num2str(jj)];
            
            if ~exist(spikeVar,'var')
                continue;
            end
            
            eval(['spikeTrain = ' spikeVar ';']);
            spikeTimes = find(spikeTrain ~= 0);
            
            if isempty(spikeTimes)
                continue;
            end
            
            spikeIndices = (spikes(:,1) == (ii-1)*64+jj-1) & (ismember(spikes(:,2),spikeTimes));
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
        end
    end
    
    writeNexFile(nexFile,outputFile);
end