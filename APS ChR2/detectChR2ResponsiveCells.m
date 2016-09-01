function detectChR2ResponsiveCells(detectionRecording,blackCondition,whiteCondition,framesPerStimulus,nTestBins,stimulusFile,saveFileSuffix,useInterspersedFlashes,useGreyFlashes,noPlot,nBlocks,repsMultiplier)
    load('channelNames.mat');
    load('EventNo.mat');
    load('spiketimestamps.mat');

    %%

%     useGreyFlashes = false;
%     useInterspersedFlashes = true;
%     detectionRecording = 2;
%     noPlot = true;
%     nTestBins = 5;

    if useInterspersedFlashes
%         blackCondition = 0;
%         whiteCondition = 49;
%         framesPerStimulus = 4;
        bw = 0.05;
%         saveFileSuffix = 'contrast_ctrl';
    else
        bw = 0.1;
    end

%     stimulusFile = 'V:\retina\John B\phd backup\aps stimuli\frontiers paper\contrast gratings\sequence.mat';
    % stimulusFile = 'V:\retina\John B\phd backup\aps stimuli\ChR2 stim test\full fields\sequence.mat';
    if exist(stimulusFile,'file')
        load(stimulusFile);
    end
    
    if nargin > 11
        nReps = nReps*repsMultiplier; %#ok<NODEF>
    else
        repsMultiplier = 1;
    end
    
    if nargin < 11
        blocks = 1:nReps;
        nBlocks = nReps;
    elseif isscalar(nBlocks)
        blocks = 1:nBlocks;
    else
        blocks = nBlocks;
        nBlocks = max(blocks);
    end
    
    if exist('nReps','var')
        assert(nBlocks <= nReps);
    end

    %%

    isBasicChR2Experiment = false;
    if useInterspersedFlashes
        fullFieldRecordings = detectionRecording;

        frameTimes = EventNo{detectionRecording}; %#ok<USENS>
        
        if exist('conditions','var')
            if isvector(conditions) %#ok<NODEF>
                conditionIndices = conditions(:);
            else
                conditionIndices = conditions(:,1);
            end
            
            nConditions = numel(conditionIndices);
            
            nStimuli = nBlocks*(nConditions/repsMultiplier+sum(~ismember([blackCondition whiteCondition],conditionIndices)));
        else
            nStimuli = size(conditionOrder,1);
        end
        
        assert(numel(frameTimes) >= nStimuli*framesPerStimulus);

        onsets = frameTimes((find(conditionOrder(1:nStimuli,1) == whiteCondition)-1)*framesPerStimulus+1);
        onsets = onsets(blocks);
        offsets = frameTimes((find(conditionOrder(1:nStimuli,1) == blackCondition)-1)*framesPerStimulus+1);
        offsets = offsets(blocks);

        nTrials = numel(onsets);
        assert(numel(offsets) == nTrials);

        tmax = 1;

        saveFile = ['ChR2_cells_' saveFileSuffix '.mat'];
    elseif useGreyFlashes
        fullFieldRecordings = detectionRecording;
        assert(numel(EventNo{detectionRecording}) == 3e4);

        onsets = EventNo{detectionRecording}(1:60:end);
        offsets = EventNo{detectionRecording}(31:60:end);

        included = conditionOrder > numel(lums)*numel(durs);

        onsets = onsets(included);
        offsets = offsets(included);

        nTrials = sum(included);

        tmax = 1;

        saveFile = 'ChR2_cells_grey.mat';
    else
        isBasicChR2Experiment = true;
        fullFieldRecordings = detectionRecording;
        detectionRecording = fullFieldRecordings(end);
        assert(numel(EventNo{detectionRecording}) == 60);

        onsets = EventNo{detectionRecording}(1:2:end-1);
        offsets = EventNo{detectionRecording}(2:2:end);

        nTrials = numel(onsets);

        tmax = 2;

        saveFile = 'ChR2_cells_ff2.mat';
    end

    nRecs = numel(fullFieldRecordings);

    if ~isBasicChR2Experiment
        assert(isscalar(detectionRecording),'Only one recording may be used for detection.');
        assert(ismember(detectionRecording,fullFieldRecordings),'Only full-fields may be used for detection.');
    end

    detectionIndex = find(fullFieldRecordings == detectionRecording);
    assert(all(fullFieldRecordings > 0 & fullFieldRecordings <= numel(EventNo)),'Missing stimulus timing file(s).');

    %%

    for ii = 1:nRecs
        mkdir(sprintf('Stim %d',fullFieldRecordings(ii)));
    end

    %%

    nSpikes = zeros(numel(spiketimestamps),1); %#ok<USENS>

    for ii = 1:nTrials
        tic;
        if useInterspersedFlashes
            nSpikesFun = @(t) sum(t > onsets(ii) & t <= onsets(ii)+tmax);
        else
            nSpikesFun = @(t) sum(t > onsets(ii) & t <= offsets(ii));
        end
    
        nSpikes = nSpikes + cellfun(nSpikesFun,spiketimestamps)';
        toc;
    end

    %%

    goodUnits = find(vertcat(channelNames{6,:}) == 0 & nSpikes >= nTrials); %#ok<USENS>
    n = numel(goodUnits);
    disp(n)

    loadNSpikesss = exist(saveFile,'file') && noPlot;
    if loadNSpikesss;
        load(saveFile,'nSpikesss','samples');
    else
        nSpikesss = cell(n,nRecs);
        samples = cell(n,nRecs);
    end
    
    p = zeros(n,nRecs,2);
    maxFr = zeros(n,nRecs,2);

    %%

    fin = fopen('originalSplit.txt');
    closeFile = onCleanup(@() fclose(fin));
    nSamples = textscan(fin, '%*s %d', 'CommentStyle', '#', 'Delimiter', ',');
    sampleRate = 7055.258405732180;
    endTimes = cumsum(double(nSamples{1}))/sampleRate;

    %%

    fig = figure;
    for ii = 1:n
        tic;
        if ~loadNSpikesss
            [fig,~,~,~,~,~,nSpikesss{ii,detectionIndex}] = rasterPlot(fig,spiketimestamps{goodUnits(ii)},[onsets; offsets],kron([1;2],ones(nTrials,1)),0,tmax,true,bw);
            unitName = channelNames{1,goodUnits(ii)};
            figFile = sprintf('Stim %d\\ffraster_channel_%s_cluster_%s',detectionRecording,unitName(3:7),unitName(8));
            saveas(fig,figFile,'fig');
            saveas(fig,figFile,'png');
        end
        
        if isBasicChR2Experiment
            if detectionRecording == 1
                t0 = 0;
            else
                t0 = endTimes(detectionRecording-1);
            end
            
            t1 = EventNo{detectionRecording}(1);
            
            spontSpikes = spiketimestamps{goodUnits(ii)};
            spontSpikes = spontSpikes(spontSpikes > t0 & spontSpikes <= t1);
            spontHist = histc(spontSpikes,t0:bw:t1);
            
            if true %~loadNSpikesss
                samples{ii,detectionIndex} = bootstrap(spontHist(1:(nTrials/2)*floor(numel(spontHist)/(nTrials/2))),nTrials/2);
            end
            
            for jj = 1:2
                maxFr(ii,detectionIndex,jj) = max(mean(nSpikesss{ii,detectionIndex}{jj}(1:2:nTrials,1:nTestBins)));
                p(ii,detectionIndex,jj) = sum(samples{ii,detectionIndex} >= maxFr(ii,detectionIndex))/numel(samples{ii,detectionIndex});
            end
        else
            if ~loadNSpikesss
                samples{ii,detectionIndex} = bootstrap(nSpikesss{ii,detectionIndex}{2},nTrials);
            end
            
            maxFr(ii,detectionIndex) = max(mean(nSpikesss{ii,detectionIndex}{1}(1:nTrials,1:nTestBins)));
            p(ii,detectionIndex) = sum(samples{ii,detectionIndex} >= maxFr(ii,detectionIndex))/numel(samples{ii,detectionIndex});
        end
        
        toc;
    end

    %%

    betterUnits = find(cellfun(@(c) sum(sum(c{1}(1:(1+isBasicChR2Experiment):nTrials,1:5))),nSpikesss(:,detectionIndex)) >= nTrials/(1+isBasicChR2Experiment));
    h = p(betterUnits,detectionIndex,1) <= 0.05/numel(betterUnits);
    bestUnits = goodUnits(betterUnits(h));
    nCells = numel(bestUnits);
    disp(nCells);

    allBestUnits = cell(nRecs,2);
    allBestUnits{detectionIndex,1} = bestUnits;
    
    betterUnits = find(cellfun(@(c) sum(sum(c{2}(1:(1+isBasicChR2Experiment):nTrials,1:5))),nSpikesss(:,detectionIndex)) >= nTrials/(1+isBasicChR2Experiment));
    h = p(betterUnits,detectionIndex,2) <= 0.05/numel(betterUnits);
    allBestUnits{detectionIndex,2} = goodUnits(betterUnits(h));

    allNSpikes = cell(nCells,nRecs);
    allNSpikes(:,detectionIndex) = nSpikesss(ismember(goodUnits,bestUnits),detectionIndex)';

    %%
    
    for ii = 1:n
        tic;
        for jj = setdiff(1:nRecs,detectionIndex)
            if ~loadNSpikesss
                [fig,~,~,~,~,~,nSpikesss{ii,jj}] = rasterPlot(fig,spiketimestamps{goodUnits(ii)},EventNo{fullFieldRecordings(jj)},repmat([1;2],30,1),0,2,true,0.1);
                unitName = channelNames{1,goodUnits(ii)};
                figFile = sprintf('Stim %d\\raster_channel_%s_cluster_%s',fullFieldRecordings(jj),unitName(3:7),unitName(8));
                saveas(fig,figFile,'fig');
                saveas(fig,figFile,'png');
            end
            
            if isBasicChR2Experiment
                if fullFieldRecordings(jj) == 1
                    t0 = 0;
                else
                    t0 = endTimes(fullFieldRecordings(jj)-1);
                end

                t1 = EventNo{fullFieldRecordings(jj)}(1);

                spontSpikes = spiketimestamps{goodUnits(ii)};
                spontSpikes = spontSpikes(spontSpikes > t0 & spontSpikes <= t1);
                spontHist = histc(spontSpikes,t0:bw:t1);
                
                if true %~loadNSpikesss
                    samples{ii,jj} = bootstrap(spontHist(1:(nTrials/2)*floor(numel(spontHist)/(nTrials/2))),nTrials/2);
                end
                
                for kk = 1:2
                    maxFr(ii,jj,kk) = max(mean(nSpikesss{ii,jj}{kk}(1:2:nTrials,1:nTestBins)));
                    p(ii,jj,kk) = sum(samples{ii,jj} >= maxFr(ii,jj,kk))/numel(samples{ii,jj});
                end
            else
                if ~loadNSpikesss
                    samples{ii,jj} = bootstrap(nSpikesss{ii,jj}{2},nTrials);
                end
                
                maxFr(ii,jj) = max(mean(nSpikesss{ii,jj}{1}(1:nTrials,1:nTestBins)));
                p(ii,jj) = sum(samples{ii,jj} >= maxFr(ii,jj))/numel(samples{ii,jj});
            end
        end
        toc;
    end
    
    %%
    
    for ii = setdiff(1:nRecs,detectionIndex)    
        for jj = 1:2
            betterUnits = find(cellfun(@(c) sum(sum(c{jj}(1:(1+isBasicChR2Experiment):nTrials,1:5))),nSpikesss(:,ii)) >= nTrials/(1+isBasicChR2Experiment));
            h = p(betterUnits,ii,jj) <= 0.05/numel(betterUnits);
            allBestUnits{ii,jj} = goodUnits(betterUnits(h));
            nCells = numel(allBestUnits{ii,jj});
            disp(nCells);
        end
        allNSpikes(:,ii) = nSpikesss(ismember(goodUnits,bestUnits),ii)';
    end

    %%

    peakFR = zeros(nCells,nRecs,2);
    allPeakFR = zeros(n,nRecs,2);

    for ii = 1:n
        tic;
        for jj = 1:nRecs
            for kk = 1:2
                allPeakFR(ii,jj,kk) = max(mean(nSpikesss{ii,jj}{kk})/0.1);
                
                idx = find(ismember(bestUnits,goodUnits(ii)));
                
                if ~isempty(idx)
                    peakFR(idx,jj,kk) = allPeakFR(ii,jj,kk);
                end
            end
        end
        toc;
    end

    %%

    save(saveFile,'-v7.3','goodUnits','betterUnits','bestUnits','allBestUnits','allNSpikes','allPeakFR','fullFieldRecordings','detectionRecording','p','maxFr','peakFR','samples','blocks','nBlocks','nSpikesss');

    %%

    close all;

    %%

    if true
        return;
    end

    %% 

    xaxislabel = 'KCl Concentration';
    xticklabels = [3 3 6 9];
    assert(numel(xticklabels) == nRecs);

    fillxs = {[2 4]};
    filllegends = {'MFA'};
    assert(numel(filllegends) == numel(fillxs));

    fillColours = [0.75 1 0.75; 0.75 0.75 1; 1 0.75 0.75; 1 1 0.75; 0.75 1 1; 1 0.75 1];

    %%

    figure;

    titlePrefixes = {'On' 'Off'};

    for ii = 1:2
        subplot(1,2,ii)
        hold on;
        uniths = plot(peakFR(:,:,ii)','Color',[0.75 0.75 0.75],'Marker','o');
        avghs = plot([mean(peakFR(:,:,ii)); median(peakFR(:,:,ii))]','LineWidth',2,'Marker','o');
        set(gca,'XTick',1:nRecs,'XTickLabel',xticklabels);
        title([titlePrefixes{ii} ' Responses']);
        xlabel(xaxislabel);

        legendhs = [uniths(1); avghs];
        legends = {'Single Units' 'Mean' 'Median'};
        yy = ylim;

        if ~isempty(fillxs)
            for jj = 1:numel(fillxs)
               h = fill(fillxs{jj}([1 2 2 1])+[-1 1 1 -1]/2,yy([1 1 2 2]),fillColours(jj,:),'EdgeColor','none');
               legendhs(end+1) = h; %#ok<AGROW>
               legends{end+1} = filllegends{jj}; %#ok<AGROW>
            end

            set(gca,'Children',[avghs; uniths; legendhs(end:-1:4)]);
        end

        xlim([1 nRecs]);
    end

    legend(legendhs,legends);
    suptitle('Peak Firing Rate');

    %%

    [~,figPrefix] = fileparts(pwd);
    figFile = sprintf([figPrefix ' peak firing rate']);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);
end