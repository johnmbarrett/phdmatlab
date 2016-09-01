function regenerateStimulus(stimulusFile,outputFile,skipStimuli)
    load(stimulusFile,'getExtraParams','getPixels','seed','textureRect','version');
    
    if ~exist('textureRect','var')
        error(['Texture dimensions unknown.  Stimulus reconstruction will be impossible for ' spikeFile]);
    end
    extraParams = validateStimulusParams(getPixels,getExtraParams);

    textureWidth = diff(textureRect([1 3]));
    textureHeight = diff(textureRect([2 4]));

    stimulus = resetStimuli(version,seed,skipStimuli);
    
    writer = VideoWriter(outputFile,'Archival');
    % TODO : frame rate
    open(writer);
    
    ii = 1;
    while ii < Inf
        tic;
        
        [pixels,extraParams,stimulus,ii] = getNextStimulus(ii,stimulus,getPixels,getExtraParams,extraParams,textureWidth,textureHeight);
        
        if any(isnan(pixels))
            break;
        end
        
        writeVideo(writer,uint8(pixels(1:extraParams.pixelSize:end,1:extraParams.pixelSize:end)));
        
        ii = ii + 1;
        
        toc;
    end
    
    close(writer);
end