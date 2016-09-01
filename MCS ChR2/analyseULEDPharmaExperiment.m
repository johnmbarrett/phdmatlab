function analyseULEDPharmaExperiment(retina,detectionRecordings)
    if nargin < 1
        retina = 'A';
    end
    
    if isnumeric(retina)
        retina = char('A'+retina-1);
    end
    
    if nargin < 2
        detectionRecordings = [3 11];
    end
    
    batchFiles = cell(16,1);
    batchFiles(setdiff(1:16,[6 15])) = repmat({'V:\retina\John B\phd backup\csharp\uledctrl\bin\Debug\stim_file_20140508_1.bat'},14,1);
    batchFiles([6 15]) = repmat({'V:\retina\John B\phd backup\csharp\uledctrl\bin\Debug\stim_file_20140710_2.bat'},2,1);
    
    [~,nTriggerss,stimSpecs,specIndices] = parseULEDArrayBatchFiles(batchFiles);
    
    concatenateSpikeTimesAndStimulusTimes(['Retina ' retina],['Retina ' retina],nTriggerss,{1:numel(specIndices)});
    
    getStimulusTimesAndConditionsFromTriggersAndStimFiles(['Retina ' retina '_vsync_times.mat'],stimSpecs,specIndices,nTriggerss,{'r' 'r' 'm'},{'full field sequence' 'flashing squares sequence' 'moving bars sequence'},['Retina ' retina]);
    
%     detectULEDResponsiveCells(retina,find(nTriggerss == 240),true,true,true,true,true);
    
%     uledSTA(retina,'flashing squares',detectionRecordings,[6 16]);
    
    movingBarAnalysisTemp(retina,detectionRecordings,1,true,true,'mbg_allresp',true,[7 17])
end