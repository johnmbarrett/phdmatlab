function luminanceStepsAnalysis(ctrlRecording,drugRecording,ctrlCellsSuffix,drugCellsSuffix,saveFileSuffix,rastersOnly,doPlot)
    fig = NaN;
    
    if nargin < 7
        doPlot = true;
    end
    
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
        drugRecording = 4;
    end
    
    if nargin < 1
        ctrlRecording = 1;
    end
    
    recordings = [ctrlRecording; drugRecording];

    %%

    load('channelNames.mat');
    load('EventNo.mat');
    load('spiketimestamps.mat');

    %%

    ctrl = load(['ChR2_cells_lumstep_' ctrlCellsSuffix '.mat']);
    drug = load(['ChR2_cells_lumstep_' drugCellsSuffix '.mat']);
    units = {ctrl.bestUnits drug.bestUnits};
    nCells = cellfun(@numel,units);
    disp(sum(nCells));

    %%

    load('./lumstep sequence.mat');
    nLums = numel(lums);
    lumstepStimuli = find(conditionOrder > 0 & conditionOrder <= nLums);

    %%
    
    saveFile = sprintf('luminance_steps_responses%s',saveFileSuffix);
    loadFile = saveFile;
    
    if ~exist(loadFile,'file');
        loadFile = 'luminance_steps_responses.mat';
    end
    
    if ~doPlot && exist(loadFile,'file')
        load(loadFile,'allNSpikes');
    else
        allNSpikes = cell(1,2);

        for ii = 1:2
            stimulusTimes = EventNo{recordings(ii)}; %#ok<USENS>
            stimulusTimes = stimulusTimes(4*(lumstepStimuli-1)+1);

            allNSpikes{ii} = cell(nCells(ii),1);

            for jj = 1:nCells(ii)
                tic;
                [~,~,~,~,~,~,allNSpikes{ii}{jj}] = rasterPlot(fig,spiketimestamps{units{ii}(jj)},stimulusTimes,conditionOrder(lumstepStimuli),[-0.5 0 0.25 0.5],[],true,0.05,{},{'Intensity'},true); %#ok<USENS>

                if ~rastersOnly || ~doPlot % TODO : rasters and everything else?
                    continue;
                end

                figFile = sprintf('Stim %d\\lumstep_raster_channel_%s',recordings(ii),channelNames{1,units{ii}(jj)}); %#ok<USENS>
                saveas(fig,figFile,'fig');
                saveas(fig,figFile,'png');
                toc;
            end
        end
    end
    
    if rastersOnly
        return;
    end

    %%

    nSpikes = cellfun(@(A) cell2mat(cellfun(@(B) cellfun(@(C) sum(sum(C(:,11:15),1),2)/50,B),A,'UniformOutput',false)'),allNSpikes,'UniformOutput',false);

    %%

    spontFR = arrayfun(@(S) cellfun(@(A) mean(reshape(sum(reshape(A{2},[50 5 4]),2),200,1)),S.allNSpikes),[ctrl drug],'UniformOutput',false);
    spontSTD = arrayfun(@(S) cellfun(@(A) std(reshape(sum(reshape(A{2},[50 5 4]),2),200,1)),S.allNSpikes),[ctrl drug],'UniformOutput',false);
    
    %%
    
    dSpikes = cellfun(@(E,S) 4*E-repmat(S',8,1),nSpikes,spontFR,'UniformOutput',false);

    %%

    figure;
    hold on;
    colours = 'bg';
    l = 0.874*lums/lums(end); %#ok<COLND>

    for ii = 1:2
        medianErrorbar(l',dSpikes{ii},'Color',colours(ii));
    end

    xlabel('Irradiance (mW/cm^2)');
    ylabel('# of Evoked Spikes');
    
    %%
    
    figFile = sprintf('lumstep_dspikes_vs_luminance%s',saveFileSuffix);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);

    %%

    pSpikes = cellfun(@(E,S) 100*(E./repmat(S',8,1)-1),nSpikes,spontFR,'UniformOutput',false);

    %%

    figure;
    hold on;

    for ii = 1:2
        pSpikes{ii}(isnan(pSpikes{ii})) = 0;
        medianErrorbar(l',pSpikes{ii}(:,all(isfinite(pSpikes{ii}),1)),'Color',colours(ii));
    end

    xlabel('Irradiance (mW/cm^2)');
    ylabel('% change in firing');
    
    %%
    
    figFile = sprintf('lumstep_pspikes_vs_luminance%s',saveFileSuffix);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);

    %%

    snr = cellfun(@(muS,muN,sigmaN) (muS-repmat(muN',8,1))./repmat(sigmaN',8,1),nSpikes,spontFR,spontSTD,'UniformOutput',false);

    %%

    figure;
    hold on;

    for ii = 1:2
        snr{ii}(isnan(snr{ii})) = 0;
        medianErrorbar(l',snr{ii}(:,all(isfinite(snr{ii}),1)),'Color',colours(ii));
    end

    xlabel('Irradiance (mW/cm^2)');
    ylabel('Signal-to-NoiseRatio');
    
    %%
    
    figFile = sprintf('lumstep_snr_vs_luminance%s',saveFileSuffix);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);
    
    %%
    
    save(saveFile,'allNSpikes','nSpikes','dSpikes','pSpikes','snr','spontFR','spontSTD');
end