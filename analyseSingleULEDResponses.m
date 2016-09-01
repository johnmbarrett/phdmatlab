function analyseSingleULEDResponses(stimRecording,spontRecording)
    spontDir = getAnalysisOutputDir(spontRecording);
    
    load(sprintf('%s\\%s_metaInfo.mat',spontDir,spontRecording.dataFile),'channelInfo');
    load(sprintf('%s\\SNR.mat',spontDir),'spikeRates');
    
    stimDir = getAnalysisOutputDir(stimRecording);
    
    load(sprintf('%s\\%s_vsync_times.mat',stimDir,stimRecording.dataFile),'vsyncTimes','recordingStartTime');
    
    ttls = vsyncTimes-recordingStartTime;
    
    stimuli = importdata('D:\backup\phd\matlab\stimuli\stim_file_20130917_1.txt',' ');
    
    commands = stimuli.textdata;
    commandIndices = find(strcmp('r',commands));
    
    X = stimuli.data(commandIndices,1);
    Y = stimuli.data(commandIndices,2);
    T = stimuli.data(commandIndices,5);
    
    pws = unique(T);
    
    figure;
    set(gcf, 'PaperPositionMode', 'auto', 'Renderer', 'zbuffer');
    
    for channel = 12:87
        spikeFile = sprintf('%s\\times_%s_channel_%d_MCD_trimmed_spikes.mat',stimDir,stimRecording.dataFile,channel);
        
        if ~exist(spikeFile,'file')
            continue;
        end
        
        channelIndex = find(strcmp(num2str(channel),{channelInfo.label}));
        
        load(spikeFile,'cluster_class');
        
        clusters = unique(cluster_class(:,1)); %#ok<NODEF>
        
        for ii = 1:numel(clusters)
            cluster = clusters(ii);
            
            if cluster == 0
                continue;
            end
            
            spikeTimes = cluster_class(cluster_class(:,1) == cluster,2)/100;
            
            if numel(spikeRates{channelIndex}) < cluster || spikeRates{channelIndex}(cluster) == 0 %#ok<USENS>
                rate = numel(spikeTimes)/max(spikeTimes);
            else
                rate = spikeRates{channelIndex}(cluster)*1000;
            end
            
            rates = zeros(16,16,numel(pws));
            
            for jj = 1:numel(X)
                onset = ttls(jj);
                spikes = spikeTimes(spikeTimes > onset & spikeTimes <= onset+0.1);
                
                responseRate = numel(spikes)/0.1;
                rateChange = responseRate/rate - 1;
                
                x = X(jj)+1;
                y = Y(jj)+1;
                t = find(pws == T(jj));
                
                rates(y,x,t) = rates(y,x,t) + rateChange/10;
            end
            
            clf;
            [rows,cols] = subplots(numel(pws));
            
            maxRate = max(max(max(rates)));
            minRate = min(min(min(rates)));
            
            for jj = 1:numel(pws)
                subplot(rows,cols,jj);
                
                if minRate == 0
                    C = rates(:,:,jj)/maxRate;
                else
                    C = (rates(:,:,jj)-minRate)/(maxRate-minRate);
                end
                
                surf(squeeze(C));
                caxis([0 1]);
                view(2);
                shading interp;
                
                xlim([1 16]);
                ylim([1 16]);
                title(sprintf('PW = %d ms',pws(jj)));
            end
            
            filename = sprintf('%s\\%s_channel_%d_cluster_%d_single_uled_responses',stimDir,stimRecording.dataFile,channel,cluster);
            saveas(gcf,filename,'fig');
            saveas(gcf,filename,'png');
        end     
    end
end