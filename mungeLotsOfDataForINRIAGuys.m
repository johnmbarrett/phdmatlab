dirs = {'.'};

for ii = 1:numel(dirs)
    cd(dirs{ii});
    load('./initRecordings.mat');
    
    for jj = 4:5
        recording = recordings(jj);
        fileDir = getAnalysisOutputDir(recording);
        load([fileDir '\' recording.dataFile '_psrh_abs']);
        load([fileDir '\' recording.dataFile '_spikes_concat_forceclustered_no_ignorenoise_no']);
        mungeDataForINRIAGuys;
    end
end
        