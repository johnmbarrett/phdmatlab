function contrastFrequencyGratingAnalysis(confounder,ctrlRecording,drugRecording,useAllResponders,nBlocks)
% confounder = 'contrast';
% confounder = 'frequency';

    if nargin < 4
        useAllResponders = false;
    end

    %%

    load('channelNames.mat');
    load('EventNo.mat');
    load('spiketimestamps.mat');
    load([confounder ' sequence.mat']);
    
    if nargin < 5
        nBlocks = nReps;
    end

    %%

    if strcmp(confounder,'contrast')
        L1 = lums([1 1 1 5 6 7])';
        L2 = lums([2 3 4 8 9 10])';
        C = (L2-L1)./(L2+L1);
        [Csorted,xvals] = sort(C);
        xticks = (1:6)';
        xticklabels = arrayfun(@(l) sprintf('%3.1f',l),100*Csorted,'UniformOutput',false);
        xlab = 'Michelson Contrast (%)';
        xlims = [0.5 6.5];
        xscale = 'linear';
        splits = [1 3; 4 6];
        subplotTitles = arrayfun(@(c) sprintf('Contrast = %3.1f%%',c),100*C,'UniformOutput',false);
    else
        xvals = (30./(barWidths*4*2))';
        xticks = wrev(xvals);
        xticklabels = xticks;
        xlab = 'Frequency (cpd)';
        xlims = [10^-2 10^-0.4];
        xscale = 'log';
        splits = [1 6];
        subplotTitles = arrayfun(@(f) sprintf('Frequency = %0.4f cpd',f),xvals,'UniformOutput',false);
    end

    %%

    ctrl = load(['ChR2_cells_' confounder '_ctrl.mat'],'bestUnits');
    drug = load(['ChR2_cells_' confounder '_drug.mat'],'bestUnits');
    units = cell(1,2);

    if useAllResponders
        allUnits = union(ctrl.bestUnits,drug.bestUnits);

        for ii = 1:2
            units{ii} = allUnits;
        end
    else
        units{1} = ctrl.bestUnits;
        units{2} = drug.bestUnits;
    end

    nCells = cellfun(@numel,units);

    %%

%     ctrlRecording = 2;
%     drugRecording = 5;
    recordings = [ctrlRecording drugRecording];

    gratingStimuli = find(~ismember(conditionOrder,[0 size(conditions,1)+1]));
    gratingStimuli = gratingStimuli(1:nBlocks*size(conditions,1));
    gratingOnsets = 4*(gratingStimuli-1)+1;

    stimulusTimes = cellfun(@(t) t(gratingOnsets),EventNo(recordings),'UniformOutput',false);
    stimulusConditions = conditions(conditionOrder(gratingStimuli),2:3);

    %%

    gratingResponses = arrayfun(@(n) zeros(nBlocks,8,6,n),nCells,'UniformOutput',false);

    % nFigs = max(conditions(:,3))*2;
    % figs = zeros(nFigs,1);

    % for ii = 1:nFigs
    %     figs(ii) = figure('Visible','off');
    % end

    for hh = 1:2
        for ii = 1:nCells(hh)
            tic;
            [~,~,~,~,~,~,nSpikesss] = rasterPlot(NaN,spiketimestamps{units{hh}(ii)},stimulusTimes{hh},stimulusConditions(:,1),0,0.25,true,0.25,stimulusConditions(:,2),{'Phase' 'Frequency'}); %#ok<USENS>
            gratingResponses{hh}(:,:,:,ii) = cell2mat(reshape(nSpikesss,[1 8 6]));
            toc;
        %     continue;
        %     
        %     unitName = channelNames{1,commonUnits(ii)};
        %     
        %     for jj = 1:nFigs
        %         isCtrl = mod(jj,2) == 1;
        %         
        %         if isCtrl
        %             recording = ctrlRecording;
        %         else
        %             recording = drugRecording;
        %         end
        %         
        %         frequency = ceil(jj/2);
        %         
        %         figFile = sprintf('Stim %d\\freqraster_channel_%s_cluster_%s_freq_%d_drug_%d',recording,unitName(3:7),unitName(8),frequency,1-isCtrl);
        %         saveas(figs(jj),figFile,'fig');
        %         saveas(figs(jj),figFile,'png');
        %     end
        %     toc;
        end
    end

    %%

    modulationDepth = arrayfun(@(n) zeros(n,6),nCells,'UniformOutput',false);

    for ii = 1:2
        for jj = 1:nCells(ii)
            tic;
            meanResps = squeeze(mean(gratingResponses{ii}(:,:,:,jj)));
            [~,maxPhase] = max(meanResps);
            [~,minPhase] = min(meanResps);

            maxResp = meanResps(sub2ind(size(meanResps),maxPhase,1:6));
            minResp = meanResps(sub2ind(size(meanResps),minPhase,1:6));

            md = (maxResp-minResp)./(maxResp+minResp);
            md(~isfinite(md)) = 0;

            modulationDepth{ii}(jj,:) = md;
            toc;
        end
    end

    %%

    figure;
    hold on;

    Y = cell2mat(cellfun(@mean,modulationDepth','UniformOutput',false))';
    E = cell2mat(cellfun(@std,modulationDepth','UniformOutput',false))';

    for ii = 1:size(splits,1)
        idx = splits(ii,1):splits(ii,2);
        errorbar(repmat(xvals(idx),1,2),Y(idx,:),E(idx,:));
    end

    legend({'Control' 'Drug'},'Location','Best');
    set(gca,'XScale',xscale,'XTick',xticks,'XTickLabel',xticklabels);
    xlabel(xlab);
    xlim(xlims);
    ylabel('Modulation Depth');

    %%

    figFile = sprintf('mod_vs_%s',confounder);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);

    %%

    % rotfun = @(pd) mod((0:7)+pd-4,8)+1;
    % rotatedResponses = zeros(size(gratingResponses));
    % 
    % for ii = 1:nCells
    %     for jj = 1:2
    %         for kk = 1:6
    %             [~,pd] = max(mean(gratingResponses(:,:,kk,jj,ii)));
    %             rotatedResponses(:,rotfun(pd),kk,jj,ii) = gratingResponses(:,:,kk,jj,ii);
    %         end
    %     end
    % end

    %%

    % maxResponses = mean(rotatedResponses(:,4,:,:,:));
    % good = squeeze(all(all(maxResponses > 0)));
    % normResponses = rotatedResponses(:,:,:,:,good)./repmat(maxResponses(:,:,:,:,good),[50 8 1 1 1]);
    % meanResponses = squeeze(mean(normResponses));
    % colours = 'rb';
    % 
    % for kk = 1:2
    %     figure;
    % 
    %     for jj = 1:6
    %         subplot(2,3,jj);
    %         hold on;
    % 
    %     %     for kk = 1:2
    %             errorbar(mean(meanResponses(:,jj,kk,:),4),std(meanResponses(:,jj,kk,:),[],4),'Color',colours(kk));
    %     %     end
    %     end
    % end

    %% 

    mutualInformation = arrayfun(@(n) zeros(n,6),nCells,'UniformOutput',false);

    x = kron((1:8)',ones(nBlocks,1));

    for ii = 1:2
        for jj = 1:nCells(ii)
            tic;

            for kk = 1:6
                y = reshape(gratingResponses{ii}(:,:,kk,jj),8*nBlocks,1);
                mutualInformation{ii}(jj,kk) = biasCorrect(x,y);
            end

            toc;
        end
    end

    %%

    figure;
    hold on;

    Y = cell2mat(cellfun(@mean,mutualInformation','UniformOutput',false))';
    E = cell2mat(cellfun(@std,mutualInformation','UniformOutput',false))';

    for ii = 1:size(splits,1)
        idx = splits(ii,1):splits(ii,2);
        errorbar(repmat(xvals(idx),1,2),Y(idx,:),E(idx,:));
    end

    legend({'Control' 'Drug'},'Location','Best');
    set(gca,'XScale',xscale,'XTick',xticks,'XTickLabel',xticklabels);
    xlabel(xlab);
    xlim(xlims);
    ylabel('Mutual Information (nats)');

    %%

    figFile = sprintf('info_vs_%s',confounder);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);

    %%

    % pid = pairwisePID(x,permute(reshape(gratingResponses,[400 6 2 nCells]),[4 1 2 3]));

    %%

    decoderPerformance = zeros(6,2);
    perfVTrials = cell(6,2);
    perfVUnits = cell(6,2);
    subUnits = 2.^(0:(nextpow2(min(nCells))-1));

    % for ii = 1:nCells
    %     tic;
        for jj = 1:2
            for kk = 1:6
                tic;
                y = reshape(gratingResponses{jj}(:,:,kk,:),8*nBlocks,nCells(jj));
                [fc,pvt,pvu] = simpleBayesianDecoder(x,y,5:5:40,subUnits,10);
                decoderPerformance(kk,jj) = fc;
                perfVTrials{kk,jj} = pvt;
                perfVUnits{kk,jj} = pvu;
                toc;
            end
        end
    %     toc;
    % end

    %%

    figure;
    hold on;

    Y = 100*decoderPerformance;

    for ii = 1:size(splits,1)
        idx = splits(ii,1):splits(ii,2);
        hs = plot(repmat(xvals(idx),1,2),Y(idx,:));

        if ii > 1
            for jj = 1:numel(hs)
                set(get(get(hs(jj),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        end
    end

    line(xlims,[12.5 12.5],'Color','k','LineStyle','--');
    legend({'Control' 'Drug' 'Chance'},'Location','Best');
    set(gca,'XScale',xscale,'XTick',xticks,'XTickLabel',xticklabels);
    xlabel(xlab);
    xlim(xlims);
    ylabel('Decoder Performance (% correct)');

    %%

    figFile = sprintf('bayes_vs_%s',confounder);
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);

    %%

    fileInfixes = {'trials' 'units'};
    ydatas = {perfVTrials perfVUnits};
    xdatas = {repmat([(5:5:40)'; 50],1,2) [repmat(subUnits',1,2); nCells]};
    xlabels = {'# Trials' '# Units'};
    xlimits = [0 55; 0.8 2^nextpow2(max(nCells))];

    for ii = 1:2
        figure;
        set(gcf,'Position',[0 0 1600 900]);
        xdata = xdatas{ii};
        ydata = cat(2,100*cell2mat(reshape(ydatas{ii}',[1 1 2 6])),100*repmat(reshape(decoderPerformance',[1 1 2 6]),[10 1 1 1]));

        for jj = 1:6
            subplot(2,3,jj);
            errorbar(xdata,squeeze(mean(ydata(:,:,:,jj))),squeeze(std(ydata(:,:,:,jj))));
            line(xlimits(ii,:),[12.5 12.5],'Color','k','LineStyle','--');
            title(subplotTitles{jj});
            xlim(xlimits(ii,:));

            if ii == 2
                set(gca,'XScale','log','XTick',xdata(:,1));
            end

            if jj == 4
                xlabel(xlabels{ii});
                ylabel('Decoder Performance (% correct)');        
            end

            if jj == 6
                legend({'Control' 'Drug' 'Chance'},'Location','Best');
            end
        end

        figFile = sprintf('bayes_vs_%s_%s',fileInfixes{ii},confounder);
        saveas(gcf,figFile,'fig');
        saveas(gcf,figFile,'png');
        close(gcf);
    end

    %%

    save([confounder '_responses.mat'],'gratingResponses','mutualInformation','modulationDepth','decoderPerformance','perfVTrials','perfVUnits');
end