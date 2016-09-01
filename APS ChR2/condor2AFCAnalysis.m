function condor2AFCAnalysis(trial,independentVariable,decoder,ctrlRecording,ctrlCells,drugRecording,drugCells)
    trial = str2double(trial);
    ctrlCells = str2double(ctrlCells);
    ctrlRecording = str2double(ctrlRecording);
    drugCells = str2double(drugCells);
    drugRecording = str2double(drugRecording);
    
    cells = arrayfun(@(n) 1:n,[ctrlCells drugCells],'UniformOutput',false);
    
    grating2AFCAnalysis(independentVariable,ctrlRecording,drugRecording,false,25,decoder,struct('section','main','stimulusIter',trial,'conditionIter',1:2,'cellIter',{cells}));
end