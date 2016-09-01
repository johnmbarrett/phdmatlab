function reversingGratingsAnalysis(ctrlRecording,drugRecording,ctrlCellsSuffix,drugCellsSuffix,saveFileSuffix,doPlot)
    if nargin < 6
        doPlot = true;
    end
    
    if nargin < 5
        saveFileSuffix = '';
    else
        saveFileSuffix = ['_' saveFileSuffix];
    end
    
    if nargin < 4
        drugCellsSuffix = 'drug';
    end
    
    if nargin < 3
        ctrlCellsSuffix = 'ctrl';
    end
    
    if nargin < 2
        drugRecording = 6;
    end
    
    if nargin < 1
        ctrlRecording = 3;
    end
    
    recordings = [ctrlRecording; drugRecording];

    %%

    load('channelNames.mat');
    load('EventNo.mat');
    load('spiketimestamps.mat');

    %%

    ctrl = load(['ChR2_cells_revgrat_' ctrlCellsSuffix '.mat']);
    drug = load(['ChR2_cells_revgrat_' drugCellsSuffix '.mat']);
    units = {ctrl.bestUnits drug.bestUnits};
    nCells = cellfun(@numel,units);

    %%

    load('./revgrat sequence.mat');
    revgratStimuli = find(conditionOrder > 0 & conditionOrder <= nContrasts);

    %%

    saveFile = sprintf('reversing_gratings_responses%s',saveFileSuffix);
    loadFile = saveFile;
    
    if ~exist(loadFile,'file')
        loadFile = 'reversing_gratings_responses.mat';
    end
    
    if ~doPlot && exist(loadFile,'file')
        load(loadFile,'allNSpikes');
    else
        allNSpikes = cell(1,2);
        fig = figure;
        set(fig,'Visible','off');

        for ii = 1:2
            stimulusTimes = EventNo{recordings(ii)};
            stimulusTimes = stimulusTimes(12*(revgratStimuli-1)+1);

            allNSpikes{ii} = cell(nCells(ii),1);

            for jj = 1:nCells(ii)
                tic;
                [~,~,~,~,~,~,allNSpikes{ii}{jj}] = rasterPlot(fig,spiketimestamps{units{ii}(jj)},stimulusTimes,conditionOrder(revgratStimuli),[-0.5 0 1 2 2.5],[],true,0.05,{},{'Contrast'},true);
                
                if ~doPlot
                    continue;
                end
                
                figFile = sprintf('Stim %d\\revgrat_raster_unit_%s',recordings(ii),channelNames{1,units{ii}(jj)});
                saveas(fig,figFile,'fig');
                saveas(fig,figFile,'png');
                toc;
            end
        end
    end

    %%

    f1Spikes = cellfun(@(A) cell2mat(cellfun(@(B) cellfun(@(C) sum(sum(C(:,26:30),1),2),B),A,'UniformOutput',false)'),allNSpikes,'UniformOutput',false);
    f2Spikes = cellfun(@(A) cell2mat(cellfun(@(B) cellfun(@(C) sum(sum(C(:,31:35),1),2),B),A,'UniformOutput',false)'),allNSpikes,'UniformOutput',false);

    %%

    f2f1SingleCell = cellfun(@(f1,f2) f2-f1,f1Spikes,f2Spikes,'UniformOutput',false);

    for ii = 1:2
        f2f1SingleCell{ii}(isnan(f2f1SingleCell{ii})) = 0;
    end

    %%

    figure;
    hold on;
    colours = 'bg';

    for ii = 1:2
        medianErrorbar(100*C,f2f1SingleCell{ii},'Color',colours(ii));
    end

    xlabel('Michelson Contrast (%)');
    ylabel('F2/F1 Ratio');
    
    %%
    
    figFile = sprintf('revgrat_f1f2sc_vs_contrast%s',saveFileSuffix);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);

    %%

    f2f1Retina = cell2mat(cellfun(@(f2f1) median(abs(f2f1),2),f2f1SingleCell,'UniformOutput',false));
%     f2f1Retina = cell2mat(cellfun(@(f1,f2) sum(f2,2)./sum(f1,2),f1Spikes,f2Spikes,'UniformOutput',false));
    
    %%
    
    figure;
    plot(100*C,f2f1Retina);
    xlim([0 70]);
    xlabel('Michelson Contrast (%)');
    ylabel('F2/F1 Ratio');
    
    %%
    
    figFile = sprintf('revgrat_f1f2ret_vs_contrast%s',saveFileSuffix);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);
    
    %%
    
    saveFile = sprintf('reversing_gratings_responses%s',saveFileSuffix);
    save(saveFile,'f1Spikes','f2Spikes','f2f1SingleCell','f2f1Retina','allNSpikes');
end