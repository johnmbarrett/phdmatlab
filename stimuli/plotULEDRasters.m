% TODO : generalise

figure;
% set(gcf,'Position',[-1600 20 1800 1200])
set(gcf,'Position',[20 100 600 400])

fileDir = '002_Stim 1';

pws = unique(sortedPWs);
% meanRates = zeros(86,1);
% nSpikess = zeros(10,10,7,8);
% latencies = zeros(10,10,7,8);

edges = -1:0.1:20.1;
load([fileDir '\Stim 1_uled_timings.mat']);
pulseTimes = stimulusTimes-recordingStartTime;
titlesuffix = '';
recordingLength = 30*60; % TODO : or is it???

allNs = zeros(0,7);
allPs = zeros(0,7);

ii = 1;
for channel = 12:87
%     if ~ismember(channel,responsiveCells(:,1))
%         continue;
%     end
    
    spikeFile = sprintf('%s\\times_Stim 1_channel_%d_MCD_trimmed_spikes.mat',fileDir,channel);
    if ~exist(spikeFile,'file');
        continue;
    end
    
    load(spikeFile,'cluster_class');
%     cluster = responsiveCells(ii,2);
    clusters = unique(cluster_class(:,1));
    
    for clusterIndex = 1:numel(clusters)
        cluster = clusters(clusterIndex);
        times = cluster_class(cluster_class(:,1) == cluster,2)/100;
        meanFR = numel(times)/recordingLength;
        
        clf;
%         hold on;
%         line(repmat([-1 20.1]',1,7),repmat(10:10:70,2,1),'Color','k','LineStyle','--');
%         set(gca,'Position',[0.03 0.06 0.95 0.89]);

%         meanRates(ii) = numel(times)/maxTime;
        
%         for jj = 1:70
%             pw = find(pws == sortedPWs(jj));
%             
%             for kk = 1:10
%                 t0 = pulseTimes(10*(sortIndices(jj)-1)+kk);
%                 t1 = t0+0.1;
%                 trialSpikes = times(times >= t0 & times <= t1)-t0;
%                 trialNum = mod(jj-1,10)+1;
%             
%                 if isempty(trialSpikes)
%                     nSpikess(trialNum,kk,pw,ii) = 0;
%                     latencies(trialNum,kk,pw,ii) = Inf;
%                 else
%                     nSpikess(trialNum,kk,pw,ii) = numel(trialSpikes);
%                     latencies(trialNum,kk,pw,ii) = trialSpikes(1);
%                 end
%             end
%         end
%         
%         ii = ii+1;

        hs = zeros(numel(edges),7);
        p = zeros(7,1);
        ns = zeros(100,7);
        
        for trial = 1:70
            trialIndex = sortIndices(trial);
            t0 = pulseTimes(10*(trialIndex-1)+1);
            t1 = pulseTimes(10*trialIndex)+1;
            pw = find(pws == sortedPWs(trial));
            
            for pulse = 1:10
                t = pulseTimes(10*(trialIndex-1)+pulse);
                w = sortedPWs(trial)/1000;
%                 fill([t t t+w t+w]-t0,trial-[1 0 0 1],[0.75 0.75 1],'EdgeColor','none');
                
                pulseSpikes = times(times >= t & times < t+0.1);
                
                ns(10*mod(trial-1,10)+pulse,pw) = numel(pulseSpikes);
                
                if ~isempty(pulseSpikes)
                    p(pw) = p(pw)+0.01;
                end
            end
            
            trialSpikes = times(times > t0-1 & times <= t1+1)-t0;

%             h = histc(trialSpikes,edges)';
%             
%             if ~isempty(h)
%                 hs(:,pw) = hs(:,pw) + h(:)/10;
%             end
            
%             plot(trialSpikes,(trial-0.5)*ones(size(trialSpikes)),'LineStyle','none','Marker','.');
        end
        
        allPs = [allPs; p'];
        
        meanN = mean(ns);
        semN = std(ns)/10;
%         hold on;
%         line([0 pws(end)*1.05],meanFR*[0.1 0.1],'Color','k','LineStyle','--');
%         errorbar(pws,meanN,2*semN);
        
        allNs = [allNs; meanN];
        
%         hold on;
%         plot(pws,p,'LineStyle','none','Marker','o','Color','k');
%         
%         sigfun = @(b,x) 1./(1+exp(-(b(1)*x+b(2))));
%         expfun = @(b,x) 1-exp(-(b(1)*x+b(2)));
%         [betaSig,~,~,~,mseSig] = nlinfit(pws,p,sigfun,[1 0]);
%         [betaExp,~,~,~,mseExp] = nlinfit(pws,p,expfun,[1 0]);
%         
%         if mseSig < mseExp
%             plotfun = @(x) sigfun(betaSig,x);
%             invfun = @(y) -(log(1/y-1)+betaSig(2))/betaSig(1);
%         else
%             plotfun = @(x) expfun(betaExp,x);
%             invfun = @(y) -(log(1-y)+betaExp(2))/betaExp(1);
%         end
%         
%         thresh = invfun(0.5);
%         
%         fplot(plotfun,[0 100],'Color','b');
%         titlesuffix = sprintf(' - Threshold: %4.2f ms',thresh);
%         
%         if thresh > 0 && thresh <= 100
%             line([thresh thresh],[0 0.5],'Color','k','LineStyle','--');
%             line([0 thresh],[0.5 0.5],'Color','k','LineStyle','--');
%         end
        
%         for pw = 1:7
%             subplot('Position',[0.03 0.05+0.1375*(pw-ii) 0.95 0.1]);
%             hold on;
%             maxY = max(max(hs(:,pw)),0.1);
%             maxX = pws(pw)/1000;
%             line(repmat(0:2:18,2,1),repmat([0;maxY],1,10),'Color','k','LineStyle','--');
%             fill(repmat(0:2:18,4,1)+repmat([0;maxX;maxX;0],1,10),kron([0;maxY],ones(2,10)),[0.75 0.75 1],'EdgeColor','none');
%             bar(edges,hs(:,pw),'histc');
%             
%             if pw == 1
%                 xlabel('Time after first pulse of train/s');
%                 ylabel('# Spikes per bin');
%             end
%             
%             set(gca,'XTick',-1:20);
%             xlim([-1 20.1]);
%             ylim([0 maxY]);
%             
%             textHandle = text(20.25,maxY,sprintf('Pulse Width: %dms',pws(pw)));
%             set(textHandle,'Rotation',270);
%         end
%         
%         mtit(sprintf('Channel %d cluster %d',channel,cluster));
        
%         continue;
        
%         xlim([-1 20]);
%         set(gca,'YTickLabel',unique(sortedPWs));
%         set(gca,'YTick',10:10:70);
%         xlabel('Time after first pulse of train/s');
%         ylabel('Pulse width/ms');

%         xlabel('Pulse width/ms');
%         xlim([0 pws(end)*1.05]);
%         ylabel('Number of spikes per pulse');
%         ylabel('Probability of at least one action potential');
%         ylim([0 1.1]);
        
%         title(sprintf('Channel %d cluster %d%s',channel,cluster,titlesuffix));
        
% %         continue;
        
%         set(gcf, 'PaperPositionMode', 'auto')
%         filename = sprintf('%s\\Stim 1_channel_%d_cluster_%d_uled_nspikes',fileDir,channel,cluster);
%         filename = sprintf('%s\\Stim 1_channel_%d_cluster_%d_uled_spikeprob',fileDir,channel,cluster);
%         filename = sprintf('%s\\Stim 1_channel_%d_cluster_%d_uled_raster',fileDir,channel,cluster);
%         filename = sprintf('%s\\Stim 1_channel_%d_cluster_%d_uled_hist',fileDir,channel,cluster);
%         saveas(gcf,filename,'fig');
%         saveas(gcf,filename,'png');
    end
end

clf;

meanAllNs = mean(allNs);
semAllNs = std(allNs)/sqrt(size(allNs,1));

errorbar(pws,meanAllNs,2*semAllNs);
xlabel('Pulse width/ms');
ylabel('Number of spikes per pulse');
title('Average of all cells');
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,[fileDir '\Stim 1_all_cells_uled_nspikes'],'fig');
saveas(gcf,[fileDir '\Stim 1_all_cells_uled_nspikes'],'png');

return

%%

for ii = 1:8
    figure;
    
    for jj = 1:7
        subplot(2,4,jj);
        
        boxplot(latencies(:,:,jj,ii));
    end
end