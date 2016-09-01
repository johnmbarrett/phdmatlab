function extraParams = validateStimulusParams(getPixels,getExtraParams)
    if ~isstruct(getPixels) && ~isa(getPixels,'function_handle') && ~iscell(getPixels)
        error(['Unsupported stimulus type for ' recording.dataFile]);
    end
    
    if isa(getExtraParams,'function_handle')
        extraParams = getExtraParams();
    else
        extraParams = getExtraParams;
    end
end