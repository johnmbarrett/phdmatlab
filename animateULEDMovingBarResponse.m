function animateULEDMovingBarResponse(recording,stimFile,varargin)
    [filedir,filename] = getAnalysisOutputDir(recording);
    
    stimInfo = importdata(stimFile);
    directions = unique(stimInfo.textdata(1:2:end-1,2));
    periods = unique(stimInfo.data(1:2:end-1,1));
    
    load(sprintf('%s\\%s_vsync_times.mat',filedir,filename),'recordingStartTime','vsyncTimes');
    stimulusTimes = vsyncTimes-recordingStartTime;
    
    options = getopt('stimulusoffset=0',varargin{:});
    barOnsets = stimulusTimes(options.stimulusoffset+(1:15:2400));
    
    allSpikeTimes = cell(8,8);
    
    function fn(spikeTimes,~,channelLabel,varargin)
        allSpikeTimes{str2double(channelLabel(2)),str2double(channelLabel(1))} = spikeTimes;
    end

    forEachChannel(recording,[],true,@fn,false,true);
    
    % TODO : different widths, different speeds
    nSpikes = zeros(141,8,8,10,8,2);
    trials = zeros(8,2);
    for ii = 1:160
        d = find(strcmp(stimInfo.textdata{2*ii-1,2},directions));
        p = find(stimInfo.data(2*ii-1,1) == periods);
        n = trials(d,p) + 1;
        trials(d,p) = n;
        t = barOnsets(ii);

        for jj = 1:141
            interval = t+(jj-1)*0.01+[0 0.1];
            for kk = 1:8
                for ll = 1:8
                    if isempty(allSpikeTimes{kk,ll})
                        continue;
                    end

                    spikes = allSpikeTimes{kk,ll};
                    spikesInTrial = spikes >= interval(1) & spikes < interval(2);
                    nSpikes(jj,kk,ll,n,d,p) = sum(spikesInTrial);
                end
            end
        end
    end
    meanResponse = squeeze(mean(nSpikes,4));
    
    maxResponse = max(max(max(max(max(meanResponse)))));

    figure;
    set(gcf,'Position',[0 0 1200 900]);
    for gg = 1:8
        for hh = 1:2
            clf;
            
            for ii = 1:8
                for jj = 1:8
                    subplot(8,8,8*(jj-1)+ii);
                    plot(0:0.01:1.4,meanResponse(:,ii,jj,gg,hh));
                    xlim([0 1.5]);
                    ylim([0 maxResponse]);
                end
            end
            
            suptitle(sprintf('Direction %s Period %d',directions{gg},periods(hh)));
            
            figFile = sprintf('%s\\moving_bar_array_response_%s_direction_%s_period_%d',filedir,filename,directions{gg},periods(hh));
            saveas(gcf,figFile,'fig');
            saveas(gcf,figFile,'png');
        end
    end
    close(gcf);

    I = uint8(254*meanResponse/maxResponse);
    square = uint8(ones(32,32));
    cmap = colormap(jet(255));

    for ii = 1:8
        for jj = 1:2
            X = zeros(32*8,32*8,1,141);

            for kk = 1:141
                X(:,:,1,kk) = kron(squeeze(I(kk,:,:,ii,jj)),square);
            end

            gifFile = sprintf('%s\\moving_bar_array_response_%s_direction_%s_period_%d.gif',filedir,filename,directions{ii},periods(jj));
            imwrite(X,cmap,gifFile,'gif','DelayTime',0.02);
        end
    end
end