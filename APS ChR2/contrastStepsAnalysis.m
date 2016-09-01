function contrastStepsAnalysis(ctrlRecording,drugRecording,ctrlCellsSuffix,drugCellsSuffix,saveFileSuffix,rastersOnly)
    fig = NaN;
    
    if nargin < 6
        rastersOnly = false;
    elseif all(logical(rastersOnly))
        fig = figure('Visible','off');
    end
    
    if nargin >= 5 && ischar(saveFileSuffix) && ~isempty(saveFileSuffix)
        saveFileSuffix = ['_' saveFileSuffix];
    else
        saveFileSuffix = '';
    end
    
    if nargin < 4
        drugCellsSuffix = 'drug';
    end
    
    if nargin < 3
        ctrlCellsSuffix = 'ctrl';
    end
    
    if nargin < 2
        drugRecording = 5;
    end
    
    if nargin < 1
        ctrlRecording = 2;
    end
    recordings = [ctrlRecording; drugRecording];

    %%

    load('channelNames.mat');
    load('EventNo.mat');
    load('spiketimestamps.mat');

    %%

    ctrl = load(['ChR2_cells_constep_' ctrlCellsSuffix '.mat']);
    drug = load(['ChR2_cells_constep_' drugCellsSuffix '.mat']);
    units = {ctrl.bestUnits drug.bestUnits};
    nCells = cellfun(@numel,units);

    %%

    load('./constep sequence.mat');
    nLums = numel(lums);
    nContrasts = nLums^2;
    constepStimuli = find(conditionOrder > 0 & conditionOrder <= nContrasts);

    %%

    allNSpikes = cell(1,2);

    for ii = 1:2
        stimulusTimes = EventNo{recordings(ii)}; %#ok<USENS>
        stimulusTimes = stimulusTimes(12*(constepStimuli-1)+1);

        allNSpikes{ii} = cell(nCells(ii),1);

        for jj = 1:nCells(ii)
            tic;
            [~,~,~,~,~,~,allNSpikes{ii}{jj}] = rasterPlot(fig,spiketimestamps{units{ii}(jj)},stimulusTimes,conditionOrder(constepStimuli),[-0.5 0 1 2 2.5],[],true,0.05,{},{'Intensity'},true); %#ok<USENS>
            
            if ~rastersOnly % TODO : rasters and everything else?
                continue;
            end
            
            figFile = sprintf('Stim %d\\constep_raster_channel_%s',recordings(ii),channelNames{1,units{ii}(jj)}); %#ok<USENS>
            saveas(fig,figFile,'fig');
            saveas(fig,figFile,'png');
            toc;
        end
    end
    
    if rastersOnly
        return;
    end

    %%

    dSpikes = cellfun(@(A) cell2mat(cellfun(@(B) cellfun(@(C) 100*(sum(sum(C(:,31:35),1),2)./sum(sum(C(:,26:30),1),2)-1),B),A,'UniformOutput',false)'),allNSpikes,'UniformOutput',false);

    for ii = 1:2
        dSpikes{ii}(isnan(dSpikes{ii})) = 0;
    end
    
    %%
    
    saveFile = sprintf('contrast_steps_responses_250ms_before_after%s',saveFileSuffix);
    save(saveFile,'allNSpikes','dSpikes');

    %%

    figure;
    hold on;
    colours = 'bg';
    contrasts = 100*(lums(conditions(2:26,3))./lums(conditions(2:26,2))-1);
    [contrasts,sortIndices] = sort(contrasts);

    for ii = 1:2
        medianErrorbar(contrasts,dSpikes{ii}(sortIndices,:),'Color',colours(ii));
    end

    xlabel('Weber Contrast (%)');
    ylabel('% change in firing');
    
    %%
    
    figFile = sprintf('constep_dspikes_vs_contrast_250ms_before_after%s',saveFileSuffix);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);
end