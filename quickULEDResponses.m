function responsiveCells = quickULEDResponses(recording)
    extractMCDSpikes(recording);
    
    extractVSyncTimings(recording,'method','threshold','threshold',1);
    
    responsiveCells = analyseULEDSquareResponses(recording,recording,'stim_file_20140508','blankttls','yes','figs','no','overwrite','yes','functions',{'squaresnew'},'suffixes',1);
    
    mungeWaveClusToNex(recording,num2str(responsiveCells(:,1)),'responsive');
end