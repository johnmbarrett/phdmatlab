function analyseFullFieldFlashes(recording,varargin)
    [fileDir,filename] = getAnalysisOutputDir(recording);
    
    load(sprintf('%s\\%s_photodiode_timings.mat',fileDir,filename),'stimulusTimes','recordingStartTime');
    stimulusTimes = stimulusTimes-recordingStartTime;
    
%     t0 = stimulusTimes(1);
    
    t = (stimulusTimes(2):0.001:stimulusTimes(62))';
    kernel = normpdf(-0.25:0.001:0.25,0,0.05)';
    kernel = kernel/sum(kernel);
    
    T = ceil(min(diff(stimulusTimes(2:62)))*1000);
    
    load(sprintf('%s.mat',filename),'getExtraParams');
    params = getExtraParams{2}; %#ok<USENS>
    
    luminances = [getExtraParams{1}.colours{1} params.colours{:}];
    onStimuli = find(luminances(1:end-1) == 0 & luminances(2:end) == 255);
    offStimuli = find(luminances(1:end-1) == 255 & luminances(2:end) == 0);
    
%     % don't know the relative intensity of 0 and 255, so can't compute
%     % actual michelson contrast values, so just do it in an abstract sense
%     levels = zeros(size(luminances));
%     
%     for hh = 1:numel(luminances)
%         levels(hh) = find(ismember(unique(luminances),luminances(hh)));
%     end
%     
%     contrasts = abs(diff(levels))./(levels(1:end-1)+levels(2:end));
%     polarities = sign(diff(levels));
%     
%     X1 = [ones(60,1) repmat(luminances(2:end)',10,1)];
%     X2 = [ones(60,1) repmat(contrasts'.*polarities',10,1)];
    
    rasterFig = figure;
    set(rasterFig,'Position',[0 0 800 600]);
    
    sdfFig = figure;
    
%     regressFig = figure;
%     set(regressFig,'Position',[0 0 400 600],'Visible','off');
    
    onoffIndices = [];
    responsivenessIndices = [];
    channels = [];
    clusters = [];
    pon = [];
    poff = [];
    
    options = getopt('figs=no',varargin{:});
    
    function fn(spikeTimes,~,channelLabel,cluster,~,~)
        channels(end+1) = str2double(channelLabel);
        clusters(end+1) = cluster;
        
        spikes = spikeTimes(spikeTimes > t(1) & spikeTimes <= t(end));
        s = zeros(size(t));
        s(ceil((spikes-t(1))*1000)) = 1;
        s = conv(s,kernel,'same');
        
        ctrlSpikes = zeros(20,1);
        stimSpikes = zeros(10,2);    
        sdf = zeros(T,2);
        
        for ii = 2:61
            stimulus = mod(ii-2,6)+1;
            
            if ~ismember(stimulus,[onStimuli offStimuli])
                continue;
            end
            
            stimOn = stimulusTimes(ii);
%             stimOff = stimulusTimes(ii+1);
            polarity = find(stimulus == [onStimuli offStimuli]);

            trial = ceil((ii-1)/6);

            ctrlSpikes(2*trial+polarity-2) = sum(spikeTimes >= stimOn-1 & spikeTimes < stimOn);
            stimSpikes(trial,polarity) = sum(spikeTimes > stimOn & spikeTimes <= stimOn+1);
            
%             contrast = find(contrasts(stimulus) == unique(contrasts));
%             polarity = find(polarities(stimulus) == unique(polarities));
            
%             nSpikes(trial,contrast,polarity) = spikesInTrial;
            
            sdf(:,polarity) = sdf(:,polarity) + 1000*s(ceil((stimOn-t(1))*1000)+(1:T))/10;
        end
        
%         sample = bootstrap(nSpikes,10,10000);
%         phi = [phi; sum(sample >= max(max(mean(nSpikes))))/numel(sample)];
%         plo = [plo; sum(sample <= min(min(mean(nSpikes))))/numel(sample)];
        pon = [pon; ranksum(ctrlSpikes,stimSpikes(:,1))];
        poff = [poff; ranksum(ctrlSpikes,stimSpikes(:,2))];
        
%         backgroundSDF = reshape(sdf(end-249:end,:),250*6,1);
%         meanFR = mean(backgroundSDF);
%         stdFR = std(backgroundSDF);

        meanFR = 0;

%         sdf(:,1) = mean(sdf(:,[1 3]),2);
%         sdf(:,2) = mean(sdf(:,[4 6]),2);
%         sdf = sdf(1:end-250,[2 5]);
        
        [~,maxI] = max(abs(sdf(1:1000,:)-meanFR));
        
        maxD = sdf(sub2ind(size(sdf),maxI,(1:2)));
        maxT = maxI/1000;
        
        onoffIndex = -diff(maxD)/sum(maxD);
        
%         responsivenessIndex = std(maxD)/mean(maxD);
%         responsivenessIndex = (abs(diff(maxD([1 2])))+abs(diff(maxD([3 4]))))/sum(maxD);

%         if stdFR == 0
%             Z = zeros(size(sdf));
%         else
%             Z = (sdf-meanFR)./stdFR;
%         end
%         
%         responsivenessIndex = max(max(abs(Z)));

        onoffIndices = [onoffIndices; onoffIndex];
%         responsivenessIndices = [responsivenessIndices; responsivenessIndex];
        
        if strcmp(options.figs,'yes')
            rasterPlot(rasterFig,spikeTimes,stimulusTimes(2:6:56),ones(10,1),[-2 0 2 4 6 8 10 12 14]',[],true,0.1);

    %         set(gcf,'Visible','on');
            figfile = sprintf('%s\\full_field_psrh_%s_channel_%s_cluster_%d',fileDir,filename,channelLabel,cluster);
            saveas(rasterFig,figfile,'fig');
            saveas(rasterFig,figfile,'png');
            
            set(0,'CurrentFigure',sdfFig);
            clf;
            hold on;
            plot((1:T)/1000,sdf);
%             line([0 (T-250)/1000],backgroundSDF*[1 1],'Color','k','LineStyle','--');
            legend({'On' 'Off'});
%             legend(cellstr(num2str(relIntensity([2 4 1 5])')))

%             line([maxT; maxT],[repmat(backgroundSDF,1,2); maxD],'Color','k','LineStyle',':');
            line([maxT; maxT],[0 0; maxD],'Color','k','LineStyle',':');

            title(sprintf('ON/OFF = %3.2f',onoffIndex));
%             title(sprintf('ON/OFF = %3.2f, Responsiveness = %3.2f',onoffIndex,responsivenessIndex));

            figfile = sprintf('%s\\sdf_%s_channel_%s_cluster_%d',fileDir,filename,channelLabel,cluster);
            saveas(sdfFig,figfile,'fig');
            saveas(sdfFig,figfile,'png');

%             Y = reshape(nSpikes,60,1);
% 
%             set(0,'CurrentFigure',regressFig);
% 
%             subplot(2,1,1);
%             hold on;
%             plot(X1(:,2),Y,'LineStyle','none','Marker','o');
%             xlim([0 255]);
%             b = regress(Y,X1);
%             fplot(@(x) b(1)+b(2)*x,[0 255]);
% 
%             subplot(2,1,2);
%             hold on;
%             plot(X2(:,2),Y,'LineStyle','none','Marker','o');
%             xlim([-255 255]);
%             b = regress(Y,X2);
%             fplot(@(x) b(1)+b(2)*x,[-255 255]);
% 
%             figfile = sprintf('%s\\spikesvintensity_%s_channel_%s_cluster_%d',fileDir,filename,channelLabel,cluster);
%             saveas(regressFig,figfile,'fig');
%             saveas(regressFig,figfile,'png');
        end
    end

    forEachChannel(recording,[],true,@fn);
    
    alpha = 0.01; %5/(numel(pon)*2);
    responsive = pon <= alpha | poff <= alpha;
    OI = onoffIndices(responsive);
    
    cellTypes = nan(numel(channels),3);
    cellTypes(:,1) = channels(:);
    cellTypes(:,2) = clusters(:);
    
    cellTypes(responsive,3) = (OI > 0.2) - (OI < -0.2);
    
    save(sprintf('%s\\%s_cell_types.mat',fileDir,filename),'cellTypes');
    
    return;
    
%     figure;
%     scatter(onoffIndices,responsivenessIndices);
    
    % 1 is OFF, 2 is ON-OFF/NR, 3 is ON
    [polarities,centroids] = kmeans(OI,3,'start',[-1; 0; 1]);
    [~,sortIndices] = sort(centroids);
    polarities(isfinite(polarities)) = sortIndices(polarities(isfinite(polarities)));
    
    figure;
    hold on;
    
    colours = 'bgr';
    for ii = 1:3
        scatter(OI(polarities == ii),ones(sum(polarities == ii),1),'MarkerEdgeColor',colours(ii));
    end
    
%     valid = responsivenessIndices < 15 & polarities == 2;
%     
%     % 1 is NR, 2 is ON-OFF
%     [index,centroids] = kmeans(responsivenessIndices(valid),2);
%     [~,sortIndices] = sort(centroids);
%     index = sortIndices(index);
    
%     figure;
%     hold on;
%     scatter(onoffIndices(polarities == 1),responsivenessIndices(polarities == 1),'MarkerEdgeColor','r');
%     scatter(onoffIndices(polarities == 3),responsivenessIndices(polarities == 3),'MarkerEdgeColor','g');
%     
%     OI = onoffIndices(valid);
%     RI = responsivenessIndices(valid);
%     
%     scatter(OI(index == 1),RI(index == 1),'MarkerEdgeColor','b');
%     scatter(OI(index == 2),RI(index == 2),'MarkerEdgeColor','y');
%     
%     ylim([0 min(15,max(responsivenessIndices))]);
    
    cellTypes = nan(numel(channels),3);
    cellTypes(:,1) = channels(:);
    cellTypes(:,2) = clusters(:);
    
    cellTypes(responsive,3) = polarities-2;
    
%     cellTypes(polarities == 1,3) = -1;
%     cellTypes(polarities == 3,3) = 1;
%     cellTypes(valid,3) = (index-1).*(index-2)./(index-2);
    
%     [k,pks] = nirenbergBootstrap(onoffIndices(isfinite(onoffIndices)),[-1 1]);
end