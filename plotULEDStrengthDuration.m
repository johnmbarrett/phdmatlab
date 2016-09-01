function plotULEDStrengthDuration(recording,varargin)
    function plotData(meanData,semData,reference,ylab,tit,filename)
        clf;
        hold on;
        
        for hh = 1:4
            errorbar(durations,meanData(:,hh),semData(:,hh),'Color',colours(hh,:),'Marker','o')
        end

        legend(cellstr([num2str(voltages) repmat('V',4,1)]));

        xlabel('Pulse Width/s');
        xlim(durations([1 end]));
        set(gca,'XTick',durations);

        ylabel(ylab);

        line(xlim',[1;1]*reference,'Color',colours(5,:),'LineStyle','--');

        title(tit);

        saveas(gcf,filename,'fig');
        saveas(gcf,filename,'png');
    end

    colours = distinguishable_colors(5);

    options = getopt('plotall=true',varargin{:});
    
    safeLoadMCDLibrary;
    [~,file] = ns_OpenFile([recording.dataFile '.mcd']);
    [~,fileInfo] = ns_GetFileInfo(file);

    fileDir = getAnalysisOutputDir(recording);
    
    load([fileDir '\' recording.dataFile '_uled_timings'],'stimulusTimes','recordingStartTime');
    pulseTimes = stimulusTimes-recordingStartTime;
    trainStarts = [1; find(diff(pulseTimes) > 1)+1; numel(pulseTimes)+1];
    
    durOrder = repmat((1:4)',20,1);
    repOrder = repmat(kron((1:5)',ones(4,1)),4,1);
    voltOrder = kron((1:4)',ones(20,1));
    
    voltages = (3:6)';
    durations = 25*(1:4)/1000;
    
    rateChanges = zeros(0,4,4);
    
    for ii = 12:87
        spikeFile = [fileDir '\times_' recording.dataFile '_channel_' num2str(ii) '_MCD_trimmed_spikes.mat'];
        
        if ~exist(spikeFile,'file')
            continue;
        end
        
        load(spikeFile,'cluster_class');
        
        clusters = unique(cluster_class(:,1)); %#ok<NODEF>
        clusters = clusters(clusters > 0);
        
        for jj = 1:numel(clusters)
            cluster = clusters(jj);
            
            spikeTimes = cluster_class(cluster_class(:,1) == cluster,2)/100;
            
            if isempty(spikeTimes)
                continue;
            end
            
            meanRate = numel(spikeTimes)/fileInfo.TimeSpan;
            assert(meanRate > 0);
            spikeRates = zeros(5,4,4);
            spikeCounts = nan(5,60,4,4);
            spikeLatencies = nan(5,60,4,4);
            
            for kk = 1:numel(trainStarts)-1
                trainStart = trainStarts(kk);
                stimStart = pulseTimes(trainStart);
                trainEnd = trainStarts(kk+1)-1;
                stimEnd = pulseTimes(trainEnd);
                spikeRate = sum(spikeTimes > stimStart & spikeTimes <= stimEnd)/(stimEnd-stimStart);
                rep = repOrder(kk);
                dur = durOrder(kk);
                volt = voltOrder(kk);
                spikeRates(rep,dur,volt) = spikeRate;
                
                pulseStarts = trainStart:2:trainEnd-1;
                for ll = 1:numel(pulseStarts)
                    pulseStart = pulseTimes(pulseStarts(ll));
                    spikes = spikeTimes(spikeTimes > pulseStart & spikeTimes <= pulseStart + 25/1000)-pulseStart;
                    spikeCounts(rep,ll,dur,volt) = numel(spikes);
                    
                    if isempty(spikes)
                        spikeLatencies(rep,ll,dur,volt) = Inf;
                    else
                        spikeLatencies(rep,ll,dur,volt) = spikes(1);
                    end
                end
            end
            
            figure;
            set(gcf,'Position',[0 0 1600 800]);
            
            for kk = 1:4
                subplot(2,2,kk);
                hold on;
                
                for ll = 1:4
                    counts = spikeCounts(:,:,ll,kk);
                    counts = counts(:,1:find(~isnan(counts(1,:)),1,'last'));
                    meanCounts = squeeze(mean(counts));
                    semCounts = 2*squeeze(std(counts))/sqrt(size(counts,1));
                    
                    t = 0:2*durations(ll):3-durations(ll);
                    x = [t wrev(t)];
                    y = [meanCounts-semCounts wrev(meanCounts+semCounts)];
                    fill(x,y,colours(ll,:),'EdgeColor','none','FaceAlpha',0.1);
                    
                    plot(t,meanCounts,'Color',colours(ll,:));
                    
                    title(sprintf('Voltage = %dV',voltages(kk)));
                    xlabel('Pulse Start Time/s');
                    ylabel('# Spikes in 1st 25 ms');
                end
            end
            
            filename = [fileDir '\' recording.dataFile '_channel_' num2str(ii) '_cluster_' num2str(cluster) '_uled_spikesperpulse'];
            saveas(gcf,filename,'fig');
            saveas(gcf,filename,'png');
            close(gcf);
            
            meanRates = squeeze(mean(spikeRates));
            
            rateChanges(end+1,:,:) = 100*(meanRates-meanRate)/meanRate; %#ok<AGROW>
            
            if all(logical(options.plotall))
                figure;
                set(gcf,'Position',[0 0 800 600]);
                filename = [fileDir '\' recording.dataFile '_channel_' num2str(ii) '_cluster_' num2str(cluster) '_uled_sd'];
                plotData(meanRates,nan(size(meanRates)),meanRate,'Firing Rate/Hz',sprintf('Channel %d cluster %d',ii,cluster),filename);
                close(gcf);
            end
        end
    end
    
    figure;
    set(gcf,'Position',[0 0 800 600]);
    
    meanRateChanges = squeeze(mean(rateChanges));
    semRateChanges = squeeze(std(rateChanges))/sqrt(size(rateChanges,1));
    
    plotData(meanRateChanges,2*semRateChanges,0,'% FR Change','',[fileDir '\' recording.dataFile '_uled_sd']);
    
    close(gcf);
    
    figure;
    set(gcf,'Position',[0 0 800 800]);
    
    for ii = 1:4
        subplot(2,2,ii);
        y = squeeze(rateChanges(:,:,ii));
        x = repmat(durations,size(y,1),1);
        plot(x,y,'LineStyle','none','Marker','.');
        hold on;
        xlim(durations([1 end])+[-0.01 0.01]);
        line(xlim',[0;0],'LineStyle','--','Color','k');
        xlabel('Pulse width/s');
        ylabel('% FR Change');
        title(sprintf('Voltage = %dV',voltages(ii)));
    end
    
    filename = [fileDir '\' recording.dataFile '_uled_sd_clusters'];
    saveas(gcf,filename,'fig');
    saveas(gcf,filename,'png');
end