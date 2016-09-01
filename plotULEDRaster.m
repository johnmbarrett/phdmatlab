function plotULEDRaster(recording,varargin)
    fileDir = getAnalysisOutputDir(recording);
    
    load([fileDir '\' recording.dataFile '_uled_timings']);
    pulseTimes = stimulusTimes-recordingStartTime;
    trainStarts = [1; find(diff(pulseTimes) > 1)+1; numel(pulseTimes)+1];
    
    durOrder = repmat((1:4)',20,1);
    repOrder = repmat(kron((1:5)',ones(4,1)),4,1);
    voltOrder = kron((1:4)',ones(20,1));
    
    voltages = (3:6)';
    durations = 25*(1:4)/100;
    
    figure;
    set(gcf,'Position',[0 0 1600 400]);
    
%     colours = distinguishable_colors(5);
    
    for ii = 36 %12:87
        spikeFile = [fileDir '\times_' recording.dataFile '_channel_' num2str(ii) '_MCD_trimmed_spikes.mat'];
        
        if ~exist(spikeFile,'file')
            continue;
        end
        
        load(spikeFile,'cluster_class');
        
        clusters = unique(cluster_class(:,1)); %#ok<NODEF>
        clusters = clusters(clusters > 0);
        
        for jj = 1 %:numel(clusters)
            cluster = clusters(jj);
            
            spikeTimes = cluster_class(cluster_class(:,1) == cluster,2)/100;
            
            if isempty(spikeTimes)
                continue;
            end
            
            for kk = 1:numel(trainStarts)-1
                stimStart = pulseTimes(trainStarts(kk));
                stimEnd = pulseTimes(trainStarts(kk+1)-1);
                
                trial = 5*(durOrder(kk)-1)+repOrder(kk);
                
                subplot(1,4,voltOrder(kk));
                hold on;
                
                trainStart = pulseTimes(trainStarts(kk));
                for ll = trainStarts(kk):2:trainStarts(kk+1)-2
                    pulseStart = pulseTimes(ll)-trainStart;
                    pulseEnd = pulseTimes(ll+1)-trainStart;
                    fill([pulseStart pulseEnd pulseEnd pulseStart],trial-[1 1 0 0],[0.75 0.75 1],'EdgeColor','none');
                end
                
                raster = spikeTimes(spikeTimes > stimStart-0.5 & spikeTimes < stimEnd+0.5)-stimStart;
                
                if ~isempty(raster)
                    plot(raster,(trial-0.5)*ones(size(raster)),'LineStyle','none','Marker','.');
                end
            end
        end
    end
    
    close(gcf);
end