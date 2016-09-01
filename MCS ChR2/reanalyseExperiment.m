function reanalyseExperiment(retina,drug,concs,volts,detectionRecordings,recordingOffset,snrRecordings,barRecordings)
    if nargin < 1
        retina = 'A';
    end
    
    if isnumeric(retina)
        retina = char('A'+retina-1);
    end
    
    if nargin < 2
        drug = '18BGA';
    end
    
    if nargin < 3
        concs = [0 0 0 0 0 0 10 20 40 80 80 80 80 80 80 0]';
    end
    
    if nargin < 4
        volts = [6 7 8 9 10 8 8 8 8 6 7 8 9 10 8 8]';
    end
    
    if nargin < 5
        detectionRecordings = [3 11];
    end
    
    if nargin < 6
        recordingOffset = 0;
    end
    
    if nargin < 8
        nRecordings = 16;
        snrRecordings = setdiff(1:16,[6 15]);
        barRecordings = [6 15];
    else
        nRecordings = max(union(snrRecordings,barRecordings));
        assert(nRecordings == numel(concs));
        assert(nRecordings == numel(volts));
    end
    
    stims = cell(nRecordings,1);
    
    for ii = 1:nRecordings
        stims{ii} = sprintf('Stim %d %d %s %dV',ii+recordingOffset,concs(ii),drug,volts(ii));
    end
        
    batchFiles = cell(nRecordings,1);
    batchFiles(snrRecordings) = repmat({'V:\retina\John B\phd backup\csharp\uledctrl\bin\Debug\stim_file_20140508_1.bat'},numel(snrRecordings),1);
    batchFiles(barRecordings) = repmat({'V:\retina\John B\phd backup\csharp\uledctrl\bin\Debug\stim_file_20140710_2.bat'},numel(barRecordings),1);
    
    [~,nTriggerss,stimSpecs,specIndices,stimFileIndices] = parseULEDArrayBatchFiles(batchFiles);
    
    concatenateSpikeTimesAndStimulusTimes(stims,['Retina ' retina],nTriggerss,stimFileIndices);
    
    getStimulusTimesAndConditionsFromTriggersAndStimFiles(['Retina ' retina '_vsync_times.mat'],stimSpecs,specIndices,nTriggerss,{'r' 'r' 'm'},{'full field sequence' 'flashing squares sequence' 'moving bars sequence'},['Retina ' retina]);
    
    detectULEDResponsiveCells(retina,find(nTriggerss == 240),true,true,true,true,true);
    
    uledSTA(retina,'flashing squares',detectionRecordings,barRecordings+[0 1]);
    
    movingBarAnalysisTemp(retina,detectionRecordings,1,true,false,'mbg_allresp',true,barRecordings+[1 2])
end