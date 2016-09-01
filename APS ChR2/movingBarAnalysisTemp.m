function movingBarAnalysisTemp(retina,detectionRecordings,responderSelection,recomputeLatencies,replotRasters,fileSuffix,recomputeBayesianDecoderPerformance,barRecordings)
    if nargin < 7
        recomputeBayesianDecoderPerformance = true;
    end
    
    if nargin < 6
        fileSuffix = '';
    elseif ~isempty(fileSuffix)
        fileSuffix = ['_' fileSuffix];
    end
        
    if nargin < 5
        replotRasters = false;
    end
    
    if nargin < 4
        recomputeLatencies = true;
    end
    
    if nargin < 3
        responderSelection = 0;
    end
    
    if nargin < 2
        detectionRecordings = [3 11];
    end
    
    if nargin < 1 || isempty(retina) || all(isnan(retina))
        isMCS = false;
        filePrefix = '';
        folderPrefix = '';
    else
        isMCS = true;
        
        if isnumeric(retina)
            retina = char('A'+retina-1);
        end
        
        filePrefix = sprintf('Retina %c_',retina);
        folderPrefix = sprintf('Retina %c\\',retina);
    end
    %%

    load([filePrefix 'channelNames.mat']);
    load([filePrefix 'spiketimestamps']);
    
    if isMCS
        load('moving bars sequence.mat');
        load(['Retina ' retina '_stimulusTimes.mat']);

        EventNo = onsetss; %#ok<USENS>

        valid = cells(:,2) > 0; %#ok<NODEF>
        spiketimestamps = spiketimestamps(valid); %#ok<NODEF>
        channelNames = channelNames(:,valid); %#ok<NODEF>

        ISIframes = 1;
    else
        load('bars sequence.mat');
        load('EventNo.mat');
    end

    %%

    if nargin < 8
        ctrl_recording = 4+3*isMCS;
        drug_recording = 7+10*isMCS;
    else
        ctrl_recording = barRecordings(1);
        drug_recording = barRecordings(2);
    end
    
    recordings = [ctrl_recording drug_recording];
    
    if ~all(isfinite(EventNo{ctrl_recording})) || ~all(isfinite(EventNo{drug_recording}))
        warning('Stimulus times missing or corrupt, aborting...\n');
        return;
    end

    for ii = 1:2
        mkdir(sprintf('%sStim %d',folderPrefix,recordings(ii)));
    end

    %%

    units = cell(1,2);

    if isMCS
        load(['Retina ' retina '_uled_responsive_cells.mat']);
        
        switch responderSelection
            case 1
                units{1} = union(responsive{detectionRecordings}); %#ok<USENS>
                units{2} = units{1};
            case 0
                for ii = 1:2
                    units{ii} = responsive{detectionRecordings(ii)};
                end
            case -1
                units{1} = intersect(responsive{detectionRecordings});
                units{2} = units{1};
        end
    else
        ctrl = load('ChR2_cells_bars_ctrl.mat','bestUnits');
        drug = load('ChR2_cells_bars_drug.mat','bestUnits');

        switch responderSelection
            case 1
                allUnits = union(ctrl.bestUnits,drug.bestUnits);

                for ii = 1:2
                    units{ii} = allUnits;
                end
            case 0
                units{1} = ctrl.bestUnits;
                units{2} = drug.bestUnits;
            case -1
                allUnits = intersect(ctrl.bestUnits,drug.bestUnits);

                for ii = 1:2
                    units{ii} = allUnits;
                end
        end
    end

    nCells = cellfun(@numel,units);

    %%

    saveFile = [filePrefix 'bar_responses' fileSuffix '.mat'];
    oldSaveFile = [filePrefix 'bar_responses.mat'];
    
    if recomputeLatencies
        ifrs = arrayfun(@(n) cell(n,1),nCells,'UniformOutput',false);
        kernel = normpdf(-125:125,0,25);
        t0 = zeros(1,2);
        idxfun = cell(1,2);

        for ii = 1:2
            assert(numel(EventNo{recordings(ii)}) == ISIframes*size(conditionOrder,1));

            t0(ii) = EventNo{recordings(ii)}(1);
            idxfun{ii} = @(t) floor(1000*(t-t0(ii)))+1;
            t1 = EventNo{recordings(ii)}(end)+isMCS*2;
            T = ceil(1000*(t1-t0(ii)))+1;

            for jj = 1:nCells(ii)
                tic;
                ifr = zeros(T,1);
                ts = idxfun{ii}(spiketimestamps{units{ii}(jj)});
                ts = ts(ts > 0 & ts <= T);
                [ut,~,idx] = unique(ts);
                ifr(ut) = accumarray(idx,1);
                ifr = conv(ifr,kernel,'same');
                ifrs{ii}{jj} = ifr;
                toc;
            end
        end
    else
        if exist(saveFile,'file')
            load(saveFile);
        elseif exist(oldSaveFile,'file')
            load(oldSaveFile);
        else
            error('Asked to reanalyse existing moving bar latencies but no moving bar responses file found.');
        end
    end

    %%

    if isMCS
        pathTime = 15*period/1000;
        speeds = (1000/16)./(period/1000)';
    else
        pathLength = 664+barWidth;
        pathTime = pathLength./(speeds/4); %#ok<NODEF>
    end

    barStimuli = find(conditionOrder > 0 & conditionOrder <= size(conditions,1)); %#ok<NODEF>
    onsets = cell(2,2);
    directions = cell(2,1);
    directionLabels = {'WE'; 'NWSE'; 'NS'; 'NESW'; 'EW'; 'SENW'; 'SN'; 'SWNE'};

    if isMCS
        directionLabels = sort(directionLabels);
    end

    for ii = 1:2
        stimulusIndices = barStimuli(conditions(conditionOrder(barStimuli),3-isMCS) == ii);
        directions{ii} = conditions(conditionOrder(stimulusIndices),2-isMCS);

        for jj = 1:2
            if isMCS
                onsets{ii,jj} = onsetss{recordings(jj)}(stimulusIndices);
            else
                onsets{ii,jj} = EventNo{recordings(jj)}(150*(stimulusIndices-1)+1);
            end
        end
    end

    nStimuli = numel(barStimuli);

    %%

    if replotRasters
        nSpikesss = arrayfun(@(n) cell(n,2),nCells,'UniformOutput',false);

        for gg = 1:2
            for hh = 1:2
                fig = figure;

                for ii = 1:nCells(gg)
                    tic;
                    [fig,~,~,~,~,~,nSpikesss{gg}{ii,hh}] = rasterPlot(fig,spiketimestamps{units{gg}(ii)},onsets{hh,gg},directionLabels(directions{hh}),[-0.5 0 pathTime(hh)+[0 0.5]],[],true,0.1);
                    unitName = channelNames{1,units{gg}(ii)};
                    figFile = sprintf('%sStim %d\\barraster_channel_%s_cluster_%s_drug_%d_speed_%d',folderPrefix,recordings(gg),unitName(3:end-1),unitName(end),gg-1,hh);
                    saveas(fig,figFile,'fig');
                    saveas(fig,figFile,'png');
                    toc;
                end

                close(fig);
            end
        end
    end

    %%

    if recomputeLatencies
        trials = zeros(8,2,2);
        ttps = arrayfun(@(n) zeros(10,8,2,n),nCells,'UniformOutput',false);

        for hh = 1:2
            for ii = 1:nStimuli
                tic;
                stimulusIndex = barStimuli(ii);
                condition = conditions(conditionOrder(stimulusIndex),:);
                direction = condition(2-isMCS);
                speed = condition(3-isMCS);
                trials(direction,speed,hh) = trials(direction,speed,hh) + 1;

                if isMCS
                    tmin = EventNo{recordings(hh)}(stimulusIndex);
                else
                    tmin = EventNo{recordings(hh)}(150*(stimulusIndex-1)+1);
                end

                imin = idxfun{hh}(tmin);
                tmax = tmin + pathTime(speed);
                imax = idxfun{hh}(tmax);

                for jj = 1:nCells(hh)
                    ifr = ifrs{hh}{jj}(imin:imax);
                    
                    if all(ifr == 0)
                        ipeak = Inf;
                    else
                        [~,ipeak] = max(ifr);
                    end
                    
                    ttps{hh}(trials(direction,speed,hh),direction,speed,jj) = ipeak/1000;
                end
                toc;
            end
        end
    end

    %%

    % rotfun = @(pd) mod((0:7)+pd-4,8)+1;
    % rotatedTTPs = zeros(size(ttps));

    modulationDepth = arrayfun(@(n) zeros(n,2),nCells,'UniformOutput',false);

    for jj = 1:2
        for ii = 1:nCells(jj)    
            for kk = 1:2
                tic;
                meanTTP = squeeze(mean(ttps{jj}(:,:,kk,ii)));
                [~,nullDir] = max(meanTTP);
                [~,prefDir] = min(meanTTP);

                maxResp = meanTTP(nullDir);
                minResp = meanTTP(prefDir);

                md = (maxResp-minResp)./(maxResp+minResp);
                md(~isfinite(md)) = 0;

                modulationDepth{jj}(ii,kk) = md;
                toc;
            end
        end
    end

    %%

    figure;
    barwitherr(cell2mat(cellfun(@std,modulationDepth','UniformOutput',false)),cell2mat(cellfun(@mean,modulationDepth','UniformOutput',false)));
    legend(arrayfun(@(s) sprintf('%d {\\mu}m/s',s),speeds,'UniformOutput',false),'Location','Best');
    set(gca,'XTickLabel',{'Control' 'Drug'});
    xlabel('Condition');
    ylabel('Modulation Depth');

    %%

    saveas(gcf,[filePrefix 'bar_mod' fileSuffix],'fig');
    saveas(gcf,[filePrefix 'bar_mod' fileSuffix],'png');
    close(gcf);

    %%

    mutualInformation = arrayfun(@(n) zeros(n,2),nCells,'UniformOutput',false);
    allMutualInformation = zeros(2,2);

    warningState = warning('query');
    warning off; %#ok<WNOFF>

    for hh = 1:2
        for ii = 1:2
            for jj = 1:nCells(hh)
                tic;
                X = struct('M',int32(8),'N',int32(1));
                X.sites = struct('label',{channelNames(1,units{hh}(jj))},'recording_tag',{{'episodic'}},'time_scale',1,'time_resolution',1/7022,'si_unit','none','si_prefix',1);
                X.categories = struct('label',{{'1'};{'2'};{'3'};{'4'};{'5'};{'6'};{'7'};{'8'}},'P',int32(10));

                for kk = 1:8
                    X.categories(kk).trials = struct('start_time',0,'end_time',pathTime(ii),'Q',num2cell(int32(ones(10,1))),'list',num2cell(ttps{hh}(:,kk,ii,jj)));
                end

                [data,counts,categories] = binlessopen(X);
                warped = binlesswarp(data);
                embedded = binlessembed(warped);

                % TODO : fix
%                 [~,~,~,I] = binlessinfo(embedded,counts,categories,X.M,struct('stratification_strategy',0,'singleton_strategy',1));
                mutualInformation{hh}(jj,ii) = 0; %I.value;
                toc;
            end
            
%             Y = struct('M',int32(8),'N',int32(1));
%             Y.sites = struct('label',{{'all_cells'}},'recording_tag',{{'episodic'}},'time_scale',1,'time_resolution',1/7022,'si_unit','none','si_prefix',1);
%             Y.categories = struct('label',{{'1'};{'2'};{'3'};{'4'};{'5'};{'6'};{'7'};{'8'}},'P',int32(10));
% 
%             for kk = 1:8
%                 Y.categories(kk).trials = struct('start_time',0,'end_time',pathTime(ii)+0.05*ii,'Q',num2cell(int32(nCells(hh)*ones(10,1))),'list',mat2cell(sort(squeeze(ttps{hh}(:,kk,ii,:)),2),ones(10,1),nCells(hh)));
%             end
%             
%             [data,counts,categories] = binlessopen(Y);

            data = reshape(ttps{hh}(:,:,ii,:),80,nCells(hh));
            
%             warped = binlesswarp(data);

            warped = mat2cell(2*(reshape(tiedrank(data(:)),size(data))-0.5)/numel(data)-1,ones(80,1),nCells(hh));
            
%             embedded = binlessembed(warped);
            
%             embedded = binlessembedu(warped);
            
            counts = int32(nCells(hh)*ones(80,1));
            
            embedded = [double(counts) vertcat(warped{:})]; % preserve ordering of cells in embedding dimensions
            
            categories = int32(kron((0:7)',ones(10,1)));

%             [~,~,~,I] = binlessinfo(embedded,counts,categories,8,struct('singleton_strategy',1));
            
            allMutualInformation(hh,ii) = 0; %I.value;
        end
    end

    warning(warningState);

    %%

    figure;
    barwitherr(cell2mat(cellfun(@std,mutualInformation','UniformOutput',false)),cell2mat(cellfun(@mean,mutualInformation','UniformOutput',false)));
    legend(arrayfun(@(s) sprintf('%d {\\mu}m/s',s),speeds,'UniformOutput',false),'Location','Best');
    set(gca,'XTickLabel',{'Control' 'Drug'});
    xlabel('Condition');
    ylabel('Mutual Information (bits)');

    %%

    saveas(gcf,[filePrefix 'bar_info' fileSuffix],'fig');
    saveas(gcf,[filePrefix 'bar_info' fileSuffix],'png');
    close(gcf);

    %%

    figure;
    bar(allMutualInformation);
    legend([arrayfun(@(s) sprintf('%d {\\mu}m/s',s),speeds,'UniformOutput',false) {'Chance'}],'Location','Best');
    set(gca,'XTickLabel',{'Control' 'Drug'});
    xlabel('Condition');
    ylabel('Mutual Information (bits)');

    %%

    saveas(gcf,[filePrefix 'all_cells_bar_perf' fileSuffix],'fig');
    saveas(gcf,[filePrefix 'all_cells_bar_perf' fileSuffix],'png');
    close(gcf);

    %%
    
    if recomputeBayesianDecoderPerformance
        decoderPerformance = zeros(2,2);
        perfVUnits = cell(2,2);
        subUnits = 2.^(0:nextpow2(min(nCells))-1);

        x = kron((1:8)',ones(10,1));

        for jj = 2:-1:1
            for kk = 1:2
                tic;
                y = reshape(ttps{jj}(:,:,kk,:),80,nCells(jj));
                [fc,~,pvu] = simpleBayesianDecoder(x,y,[],subUnits,10,'mixedBernoulliGaussian',0,pathTime(kk)+0.05*kk);
                decoderPerformance(jj,kk) = fc;
                perfVUnits{jj,kk} = pvu;
                toc;
            end
        end

        perfVUnits = cell2mat(reshape(perfVUnits,[1 1 2 2]));

        %%

        figure;
        bar(100*decoderPerformance);
        line(xlim,[12.5 12.5],'Color','k','LineStyle','--');
        legend([arrayfun(@(s) sprintf('%d {\\mu}m/s',s),speeds,'UniformOutput',false) {'Chance'}],'Location','Best');
        set(gca,'XTickLabel',{'Control' 'Drug'});
        xlabel('Condition');
        ylabel('Decoder Performance (%)');

        %%

        saveas(gcf,[filePrefix 'bar_perf' fileSuffix],'fig');
        saveas(gcf,[filePrefix 'bar_perf' fileSuffix],'png');
        close(gcf);

        %%

        figure;
        set(gcf,'Position',[680 558 800 420]);

        xx = [0.8 2^nextpow2(max(nCells))];

        for jj = 1:2
            subplot(1,2,jj);
            errorbar([repmat(subUnits',1,2);nCells],[squeeze(mean(100*perfVUnits(:,:,:,jj))); 100*decoderPerformance(:,jj)'],[squeeze(std(100*perfVUnits(:,:,:,jj))); 0 0]);
            line(xx,[12.5 12.5],'Color','k','LineStyle','--');
            set(gca,'XScale','log','XTick',[subUnits max(nCells)]);
            title(sprintf('Speed = %d {\\mu}m/s',speeds(jj)));
            xlim(xx);
        end

        legend({'Control' 'Drug' 'Chance'},'Location','Best');

        subplot(1,2,1);
        xlabel('# Cells');
        ylabel('Decoder Performance (%)');

        %%

        saveas(gcf,[filePrefix 'bar_perf_vs_units' fileSuffix],'fig');
        saveas(gcf,[filePrefix 'bar_perf_vs_units' fileSuffix],'png');
        close(gcf);
    end

    %%

    save(saveFile,'-v7.3','ifrs','ttps','mutualInformation','allMutualInformation','modulationDepth','decoderPerformance','perfVUnits');
end