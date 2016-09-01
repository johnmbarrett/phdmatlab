function mungeWaveClusToNexForAnderson(filenames,outputFile,samplingFrequency)
    if nargin < 3
        samplingFrequency = 25000;
    end

    nexFile = nexCreateFileData(samplingFrequency);
    
    for ii = 1:numel(filenames)
        load(filenames{ii},'index','spikes')
        nexFile = nexAddWaveform(nexFile,25000,index(:)/100,1e3*spikes',filenames{ii}); %#ok<NODEF>
    end
    
    if ~strcmp(outputFile(end-3:end),'.nex')
        outputFile = [outputFile '.nex'];
    end
    
    writeNexFile(nexFile,outputFile);
end