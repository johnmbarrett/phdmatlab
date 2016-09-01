fn = @(spikeTimes,channelIndex,channelLabel,cluster,maxClusters,spikes) rasterPlot(spikeTimes,(5:4:41)',ones(10,1),[-0.5 0 1 2 2.5],[],false,NaN);

[figs,channels,clusters] = forEachChannel('Stim 9 Lucivid max gain full fields',[],true,fn,1,true);

%%

for ii = 1:numel(figs)
    filename = sprintf('Stim 9 Lucivid max gain full fields\\raster_channel_%s_cluster_%d',channels{ii},clusters(ii));
    saveas(figs{ii},filename,'fig');
    saveas(figs{ii},filename,'png');
end

% 19/10/2014 - analysing mutual info for LED moving bars

% for aa = 1:9
%     expDir = sprintf('JBOG%04d',experiments(aa));
%     cd(expDir);
%     
%     try
%         for bb = 1:2
%             ctrlIndex = 6+16*(bb-1)+stimOffsets(aa);
%             drugIndex = ctrlIndex+9;
%             
%             for cc = 1:3
%                 ctrlRec = sprintf('Stim %d 0 %s %dV',ctrlIndex,drugs{cc},middleVolts(aa,cc));
%                 drugRec = sprintf('Stim %d 80 %s %dV',drugIndex,drugs{cc},middleVolts(aa,cc));
%                 
%                 if ~(exist(ctrlRec,'dir') && exist(drugRec,'dir'))
%                     continue;
%                 end
%                 
%                 rec = {ctrlRec drugRec};
%                 
%                 for dd = 1:2
%                     psrhfile = sprintf('%s\\%s_moving_bar_psrh.mat',rec{dd},rec{dd});
%                     
%                     if ~exist(psrhfile,'file')
%                         warning('Missing PSRH file #%d for experiment %d drug %s',dd,aa,drugs{cc});
%                         continue;
%                     end
%                     
%                     load(psrhfile);
%                     
%                     for ii = 1:8
%                         for jj = 1:2
%                             for kk = 1:10
%                                 stimulusTimings{ii,jj,1}(2,1,kk) = stimulusTimings{ii,jj,1}(1,1,kk)+0.75*jj; %#ok<SAGROW>
%                             end
%                         end
%                     end
%                     
%                     save(psrhfile,'stimulusTimings','-append');
%                     
%                     concatenateSpikes(rec,'forceclustered','no','ignorenoise','yes');
%                 end
%                 
%                 respfile = sprintf('%s_%s_responsive_cells.mat',expDir,drugs{cc});
%                 
%                 if ~exist(respfile,'file')
%                     warning('Missing responsive cells file for experiment %d drug %s',experiments(aa),drugs{cc});
%                     continue;
%                 end
%                 
%                 analyseGollischMeisterResponses({ctrlRec drugRec},'latencytype','peak','significance','file','rasterfilesuffix','moving_bar_psrh','infofilesuffix',sprintf('%s_info',drugs{cc}),'infovar',1,'confvar',2,'responsefile',respfile);
%             end
%         end
%         
%         cd ..;
%     catch err
%         logMatlabError(err);
%         cd ..;
%     end
% end

% 13/06/2014 - still trying....
% 
% indices = find(diff(squareData > 0.58) == 1);
% indices = indices([1; find(diff(indices) > 40000) + 1]);
% 
% dI = round(diff(indices)/25000);
% gaps = find(dI > 2);
% 
% inserted = 0;
% for ii = 1:numel(gaps)
%     index = gaps(ii);
%     missing = dI(index)/2-1;
%     dt = diff(indices([index index+1]+inserted));
%     dt = floor(dt/(missing+1));
%     indices = [indices(1:index+inserted); indices(index+inserted)+(1:missing)'*dt; indices(index+inserted+1:end)];
%     inserted = inserted + missing;
% end

% 12/06/2014 - trying to fix messed up timing signals from artificial eye
% LEDs

% for ii = 1:65
%     indices = find(diff(barData(75000+100000*(ii-1)+(1:50000)) > thresholds(ii)) == 1);
%     indices = indices([1; find(diff(indices) > 1000) + 1]);
%     
%     if gaps(ii) > 0
%         indices = [indices(1:gaps(ii)-1); mean(indices(gaps(ii)-[1 0])); indices(gaps(ii):end)];
%     end
%     
%     meanDiff = mean(diff(indices));
%     indices = [indices(1)-(detectable(ii)-1:-1:1)'*meanDiff; indices];
%     
%     if numel(indices) < 14
%         indices = [indices; indices(end)+(1:14-numel(indices))'*meanDiff];
%     end
%     
%     barIndices(ii,:) = indices(1:14)';
% end
% 
% for ii = 66:118
%     indices = find(diff(barData(85000+100000*(ii-1)+(1:50000)) > thresholds(ii)) == 1);
%     indices = indices([1; find(diff(indices) > 1000) + 1]);
%     
%     if gaps(ii) > 0
%         indices = [indices(1:gaps(ii)-1); mean(indices(gaps(ii)-[1 0])); indices(gaps(ii):end)];
%     end
%     
%     meanDiff = mean(diff(indices));
%     indices = [indices(1)-(detectable(ii)-1:-1:1)'*meanDiff; indices];
%     
%     if numel(indices) < 14
%         indices = [indices; indices(end)+(1:14-numel(indices))'*meanDiff];
%     end
%     
%     barIndices(ii,:) = indices(1:14)';
% end
% 
% for ii = 119:160
%     indices = find(diff(barData(95000+100000*(ii-1)+(1:50000)) > thresholds(ii)) == 1);
%     indices = indices([1; find(diff(indices) > 1000) + 1]);
%     
%     if gaps(ii) > 0
%         indices = [indices(1:gaps(ii)-1); mean(indices(gaps(ii)-[1 0])); indices(gaps(ii):end)];
%     end
%     
%     meanDiff = mean(diff(indices));
%     indices = [indices(1)-(detectable(ii)-1:-1:1)'*meanDiff; indices];
%     
%     if numel(indices) < 14
%         indices = [indices; indices(end)+(1:14-numel(indices))'*meanDiff];
%     end
%     
%     barIndices(ii,:) = indices(1:14)';
% end

% % 06/06/2014 - finding cell RF from moving bar
% width = 2;
% pathLength = 16-2;
% 
% pixels = zeros(16,16,pathLength);
% 
% theta = 1*pi/4;
% 
% [Y,X] = ndgrid(0:15,0:15);
% 
% for ii = 0:pathLength-1
%     w = (X-8)*cos(theta)-(Y-8)*sin(theta)+8;
%     pixels(:,:,ii+1) = (w >= ii & w < ii + width);
% end
% 
% figure
% for ii = 1:pathLength
%     surf(X,Y,255*pixels(:,:,ii));
%     colormap(gray);
%     view(2);
%     input('...');
% end

% % 27/05/2014 - visualising array-wide responses to moving bars
% nSpikes = zeros(141,8,8,10,8,2);
% trials = zeros(8,2);
% for ii = 1:160
%     d = find(strcmp(A.textdata{2*ii-1,2},directions));
%     p = find(A.data(2*ii-1,1) == periods);
%     n = trials(d,p) + 1;
%     trials(d,p) = n;
%     t = barOnsets(ii);
%     
%     for jj = 1:141
%         interval = t+(jj-1)*0.01+[0 0.1];
%         for kk = 1:8
%             for ll = 1:8
%                 if isempty(allSpikeTimes{kk,ll})
%                     continue;
%                 end
%                 
%                 spikes = allSpikeTimes{kk,ll};
%                 spikesInTrial = spikes >= interval(1) & spikes < interval(2);
%                 nSpikes(jj,kk,ll,n,d,p) = sum(spikesInTrial);
%             end
%         end
%     end
% end
% meanResponse = squeeze(mean(nSpikes,4));
% %%
% maxResponse = max(max(max(max(max(meanResponse)))));
% 
% for gg = 1:8
%     for hh = 1:2
%         figure
%         for ii = 1:8
%             for jj = 1:8
%                 subplot(8,8,8*(jj-1)+ii);
%                 plot(0:0.01:1.4,meanResponse(:,ii,jj,gg,hh));
%                 xlim([0 1.5]);
%                 ylim([0 maxResponse]);
%             end
%         end
%         suptitle(sprintf('Direction %s Period %d',directions{gg},periods(hh)));
%     end
% end
% 
% %%
% I = uint8(254*meanResponse/maxResponse);
% square = uint8(ones(32,32));
% cmap = colormap(jet(255));
% 
% for ii = 1:8
%     for jj = 1:2
%         X = zeros(32*8,32*8,1,141);
%         
%         for kk = 1:141
%             X(:,:,1,kk) = kron(squeeze(I(kk,:,:,ii,hh))',square);
%         end
%         
%         imwrite(X,cmap,sprintf('response_d%d_p%d.gif',ii,jj),'gif','DelayTime',0.02);
%     end
% end

% fudge = [1 -1 -1 -1 1 1 1 1];
% for ii = 1:8
% x0 = startX(ii);
% y0 = startY(ii);
% theta = (ii-1)*pi/4;
% m = round(tan(theta+pi/2));
% for jj = 1:16
% x = x0+(1-2*(x0 > 0))*jj;
% y = y0+(1-2*(y0 > 0))*jj;
% I = (fudge(ii)*(Y-y-m.*(X-x)) <= 0);
% J = (fudge(ii)*(Y-(y+width*sin(theta))-m*(X-(x+width*cos(theta)))) > 0);
% surf(X,Y,255*(I&J)); view(2);
% % pause
% end
% end

% [Y,X] = ndgrid(-7:8,-7:8);
% 
% exactTrig = @(f,x) sign(f(x)).*sqrt(round(2*f(x).^2)/2);
% 
% for ii = 1:8
%     theta = (ii-1)*pi/4;
%     
%     W = X*exactTrig(@cos,theta)-Y*exactTrig(@sin,theta)+8;
% %     Z = X*sin(theta)+Y*cos(theta)+8;
%     
%     for jj = 1:16
%         I = W >= jj & W < jj+width;
%         surf(X,Y,255*I);
%         view(2);
%     end
% end