function computeFullFieldPSTH(recording,varargin)
    [fileDir,filename] = getAnalysisOutputDir(recording);
    
    stimTimesFile = sprintf('%s\\%s_photodiode_timings.mat',fileDir,filename);
    
    if ~exist(stimTimesFile,'file')
        error('Could not find stim file for recording %s\n',filename);
    end
    
    options = getopt('bw=0.01 ffindices=NaN stimlength=2',varargin{:});
    
    load(stimTimesFile,'stimulusTimes','recordingStartTime');
    
    if isnan(options.ffindices)
        options.ffindices = 2:2:40;
    end
    
    nTrials = numel(options.ffindices)/2;
    
    stimulusTimes = stimulusTimes(options.ffindices) - recordingStartTime;
    
    psth = zeros(ceil(2*options.stimlength/options.bw),0);
    cells = zeros(0,2);
    baseline = [];
    
    preStimInterval = stimulusTimes(1)-[5 0];
    
    function fn(spikeTimes,~,channelLabel,cluster,varargin)
        cells(end+1,:) = [str2double(channelLabel) cluster];
        
        baseline = [baseline; sum(spikeTimes >= preStimInterval(1) & spikeTimes < preStimInterval(2))];
        
        [fig,~,~,hists] = rasterPlot(spikeTimes,stimulusTimes,repmat({'On';'Off'},nTrials,1),0,options.stimlength,true,options.bw,[],{'Polarity'},false);
        hists = hists/nTrials;
        
        figFile = sprintf('%s\\ffpsth_%s_channel_%s_cluster_%d',fileDir,filename,channelLabel,cluster);
        saveas(fig,figFile,'fig');
        saveas(fig,figFile,'png');
        close(fig);
        
        % swap columns because when we sort the conditions Off comes before On
        psth = [psth reshape(hists(:,[2 1]),size(psth,1),1)];
    end

    forEachChannel(recording,[],true,@fn,true,true);
    
    baseline = baseline/(5/options.bw);
    
    biasIndices = zeros(size(cells,1),1);
    
    for ii = 1:size(cells,1);
        p = psth(:,ii);
        [~,~,b] = g_biasIndex([p g_gausssmooth(p,5)],floor(size(psth,1)/2),baseline(ii),13);
        biasIndices(ii) = b;
    end
    
    respIndices = find(biasIndices ~= 2);
    respCells = cells(respIndices,:); %#ok<NASGU,FNDSB>
    
    save(sprintf('%s\\%s_ffpsth',fileDir,filename),'psth','cells','baseline','respIndices','respCells');
end
