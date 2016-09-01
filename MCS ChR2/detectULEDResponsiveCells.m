function detectULEDResponsiveCells(retina,fullFieldRecordings,recomputeRasters,plotRasters,recomputeThresholds,plotThresholds,recomputeSNR)
    if nargin < 1
        retina = 'A';
    end
    
    if nargin < 2
        fullFieldRecordings = setdiff(1:18,[6 7 16 17]);
    end
    
    if nargin < 3
        recomputeRasters = true;
    end
    
    if nargin < 4
        plotRasters = true;
    end
    
    if nargin < 5
        plotThresholds = true;
    end
    
    if isnumeric(retina)
        retina = char('A'+retina-1);
    end

    %%

    sequence = load('full field sequence.mat');
    conditions = sequence.conditions;
    conditionOrder = sequence.conditionOrder;
    pws = sequence.duration;

    nStimuli = numel(conditionOrder);
    nConditions = size(conditions,1);
    nBlocks = nStimuli/nConditions;

    trainIndices = 1:2:nBlocks;
    testIndices = 2:2:nBlocks;
    nTests = numel(testIndices);

    %%

    load(['Retina ' retina '_spiketimestamps.mat']);
    spikeTimess = spiketimestamps(cells(:,2) > 0); %#ok<NODEF>
    cells = cells(cells(:,2) > 0,:);
    nCells = numel(spikeTimess);

    %%

    load(['Retina ' retina '_stimulusTimes.mat']);

    nRecordings = numel(fullFieldRecordings);

    for ii = 1:nRecordings
        % TODO : this is wrong and loses condition information
        mkdir(sprintf('Retina %c\\Stim %d',retina,fullFieldRecordings(ii)));
    end

    bws = [0.01 0.1];
    ts = [-0.1 0 0.1 0.2; -1 0 0.1 1] ;
    prefixes = {'shortraster' 'longraster'};

    %%

    allBinss = cell(nCells,nRecordings);
    spontSpikess = cell(nCells,nRecordings);
    stimSpikess = cell(nCells,nRecordings);
    meds = zeros(nCells,nRecordings);
    mus = zeros(nCells,nRecordings);
    sigmas = zeros(nCells,nRecordings);
    ps = zeros(nCells,nRecordings);
    nSpikes = zeros(nTests,nConditions,nCells,nRecordings);
    pSpikes = zeros(nConditions,nCells,nRecordings);
    snr = zeros(nConditions,nCells,nRecordings);
    responsive = cell(nRecordings,1);
    thresholds = Inf(nCells,nRecordings);

    saveFile = ['Retina ' retina '_uled_responsive_cells.mat'];

    if ~recomputeRasters || ~recomputeThresholds || ~recomputeSNR
        load(saveFile);
    end

    for ii = 1:nRecordings
        recIndex = fullFieldRecordings(ii);
        stimulusTimes = onsetss{recIndex}; %#ok<USENS>

        for jj = 1:nCells
            tic;
            channel = cells(jj,1);
            cluster = cells(jj,2);

            if recomputeRasters
                for kk = 1:2
                    [figs,~,~,~,~,~,allBins] = rasterPlot(spikeTimess{jj},stimulusTimes,pws(conditionOrder),ts(kk,:),[],true,bws(kk),ones(size(conditionOrder)),{'PW'});

                    if kk == 2
                        allBinss{jj,ii} = cat(3,allBins{:});
                    end

                    if ~plotRasters
                        close(figs);
                        continue;
                    end

                    for ll = 1:numel(figs)
                        fig = figs(ll);
                        set(fig,'Position',[0 0 1600 900]);
                        prefix = prefixes{kk};
                        filename = sprintf('Retina %c\\Stim %d\\ff%s_channel_%d_cluster_%d',retina,recIndex,prefix,channel,cluster);
                        saveas(fig,filename,'fig');
                        saveas(fig,filename,'png');
                        close(fig);
                    end
                end
            end
            
            if ~recomputeSNR
                toc;
                continue;
            end

            spontSpikes = reshape(allBinss{jj,ii}(:,1:10,:),nBlocks*nConditions*10,1);
            spontSpikess{jj,ii} = spontSpikes;
            medianSpontSpikes = median(spontSpikes);
            meds(jj,ii) = medianSpontSpikes;

            mu = mean(spontSpikes);
            mus(jj,ii) = mu;

            % 0 std => Inf SNR, so put a lower bound on std equivalent to one
            % spontaneous spike in the entire recording
    %         sigma = max(std(spontSpikes),sqrt(1/numel(spontSpikes)));

            % actually, if I'm using ranks instead of absolute values, I can
            % probably get away with including Infs
            sigma = std(spontSpikes);
            sigmas(jj,ii) = sigma;

            stimSpikes = squeeze(allBinss{jj,ii}(:,11,:));
            stimSpikess{jj,ii} = stimSpikes;

    %         if reBootstrap
            sample = bootstrap(spontSpikes,numel(trainIndices),10000,true); % or false?
            maxStimSpikes = stimSpikes(trainIndices,end);

            ps(jj,ii) = sum(sample >= median(maxStimSpikes))/numel(sample); % or mean?
    %         end

            nSpikes(:,:,jj,ii) = stimSpikes(testIndices,:);
            pSpikes(:,jj,ii) = sum(nSpikes(:,:,jj,ii) > medianSpontSpikes)/nTests; % or mu?
            snr(:,jj,ii) = mean(nSpikes(:,:,jj,ii)-mu)/sigma;
            snr(isnan(snr(:,jj,ii)),jj,ii) = 0;
            toc;
        end

        if recomputeSNR
            responsive{ii} = find(fdrcorrect(ps(:,ii),0.05));
        end
        
        if ~recomputeThresholds
            continue;
        end

        if plotThresholds
            figure;
        end

        sigmoid = @(b,x) 1./(1+exp(-(x-b(1))/b(2)));
        beta0 = [0 1];

        for jj = 1:nCells
            tic;
            % TODO : mean?
            bp = nlinfit(pws,pSpikes(:,jj,ii),sigmoid,beta0);
            thresholds(jj,ii) = max(0,bp(1));

            if ~plotThresholds
                continue;
            end

            clf;
            hold on;

            plot(pws,pSpikes(:,jj,ii),'Color','k','LineStyle','none','Marker','o');

            fplot(@(x) sigmoid(bp,x),xlim,'Color','b');

            title(sprintf('Channel %d Cluster %d - Threshold: %5.2f ms',cells(jj,1),cells(jj,2),thresholds(jj,ii)));
            xlabel('Pulse Width/ms');
            ylabel('Response Probability');
            ylim([0 1]);

            line(xlim,[0.5 0.5],'Color','k','LineStyle','--');
            line(bp([1 1]),[0 1],'Color','k','LineStyle','--');

            figFile = sprintf('Retina %c\\Stim %d\\thresholdsnew_channel_%d_cluster_%d',retina,recIndex,cells(jj,1),cells(jj,2));
            saveas(gcf,figFile,'fig');
            saveas(gcf,figFile,'png');
            toc;
        end
        
        if plotThresholds
            close(gcf);
        end
    end
    
    if recomputeThresholds || recomputeRasters || recomputeSNR
        save(saveFile,'-v7.3','allBinss','spontSpikess','stimSpikess','meds','mus','sigmas','ps','nSpikes','pSpikes','snr','responsive','thresholds','cells');
    end

    %%

    responseIndices = intersect(responsive{[3 11]});
    nResponses = numel(responseIndices);
    datas = {nSpikes nSpikes-repmat(reshape(mus,[1 1 size(mus)]),[size(nSpikes,1) size(nSpikes,2) 1 1]) reshape(snr,[1 size(snr)])};
    permutations = {[3 4 5 1 2] [3 4 5 1 2] [3 2 4 5 1]};
    exptIndicess = {[1:5; 9:13]' [3 6 7 8 11 14]' [3 11]};
    pwIndicess = {6 6 (1:6)'};
    xlabels = {'Voltage' 'Concentration' 'Flash Duration'};
    xticklabels = {6:10 {'Ctrl' '10 uM' '20 uM' '40 uM' '80 uM' 'Wash'} pws'};
    xvals = {(6:10)' (1:6)' pws};
    ylabels = {'# Spikes' '# Extra Spikes' 'SNR'};
    figSuffixes = {'volt' 'conc' 'dur'};

    %%

    for ii = 1:3
        figure;
        set(gcf,'Position',[0 0 1920 900]);
        exptIndices = exptIndicess{ii};
        pwIndices = pwIndicess{ii};

        for jj = 1:3
            subplot(1,3,jj);
            
            data = permute(reshape(mean(datas{jj}(:,pwIndices,responseIndices,exptIndices),1),[1 numel(pwIndices) nResponses size(exptIndices,1) size(exptIndices,2)]),permutations{ii});
            medData = squeeze(median(data,1));
            
            errorbar(repmat(xvals{ii},1,size(exptIndices,2)),medData,medData-squeeze(prctile(data,25,1)),squeeze(prctile(data,75,1))-medData);

            set(gca,'XTick',xvals{ii}','XTickLabel',xticklabels{ii});

            if jj == 1
                xlabel(xlabels{ii});
            end

            ylabel(ylabels{jj});
        end

        if ismember(ii,[1 3])
            legend({'Control' 'Drug'},'Location','Best');
        end

        figFile = ['Retina_' retina '_uled_resps_vs_' figSuffixes{ii}];
        saveas(gcf,figFile,'fig');
        saveas(gcf,figFile,'png');
        close(gcf);
    end
end