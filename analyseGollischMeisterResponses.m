function analyseGollischMeisterResponses(recordings,varargin)
    opts = getopt( ...
            ['latencytype=''firstspike'' ' ...
             'significance=''rf'' ' ...
             'rasterfilesuffix=''psrh_abs'' ' ...
             'infofilesuffix=''info'' ' ...
             'infovar=2 ' ...
             'confvar=1 ' ...
             'psrhfile=NaN ' ...
             'spiketimesfile=NaN ' ...
             'responsefile=NaN'
            ],varargin{:});
        
    nRecordings = numel(recordings);
    for ee = 1:nRecordings
        recording = recordings{ee};
        %% INTITIALIZATION and defining significant responses
        if ischar(opts.psrhfile)
            psrhFile = opts.psrhfile;
            
            [fileDir,filename] = fileparts(psrhFile);
            
            if isempty(fileDir)
                fileDir = '.';
            end
        else
            [fileDir,filename] = getAnalysisOutputDir(recording);
            
            psrhFile = [fileDir '\' filename '_' opts.rasterfilesuffix '.mat'];
        end
        
        load(psrhFile,'cells','levels','valuess','rasters','stimulusTimings','histograms','edgess');
        
        if ee == 1 && iscell(cells) %#ok<NODEF>
            cellIndices = struct('maxIndex',0);
        end

        permutation = [1 opts.confvar+1 opts.infovar+1 4];
        rasters = permute(rasters,permutation);
        histograms = permute(rasters,permutation); %#ok<NASGU>
        edgess = permute(rasters,permutation); %#ok<NASGU>

        stimulusTimings = permute(stimulusTimings,[opts.confvar opts.infovar 3]);

        nCells = size(cells,1);

        if iscell(cells)
            for ii = 1:nCells
                channel = cells{ii,1};
                cluster = cells{ii,2};

                if ~isfield(cellIndices,channel)
                    indices = zeros(1,cluster);
                    indices(cluster) = cellIndices.maxIndex+1;
                    cellIndices.(channel) = indices;
                elseif numel(cellIndices.(channel)) < cluster || cellIndices.(channel)(cluster) == 0
                    cellIndices.(channel)(cluster) = cellIndices.maxIndex+1;
                end
                
                cellIndices.maxIndex = cellIndices.maxIndex+1;
            end
                    
        else
            if ee == 1
                cellIndices = zeros(56,56,max(cells(:,3))+1);
            elseif max(cells(:,3)+1) > size(cellIndices,3)
                cellIndices(:,:,end+1:max(cells(:,3))+1) = 0;
            end

            index = max(max(max(cellIndices)));
            for ii = 1:nCells
                if cellIndices(cells(ii,1),cells(ii,2),cells(ii,3)+1) == 0
                    index = index+1;
                    cellIndices(cells(ii,1),cells(ii,2),cells(ii,3)+1) = index;
                end
            end
        end

        if ischar(opts.spiketimesfile) && exist(opts.spiketimesfile,'file')
            sortedSpikeTimes = load(opts.spiketimesfile,'spiketimestamps');
            sortedSpikeTimes = sortedSpikeTimes.spiketimestamps;
        else
            fileOrder = {'No' 'Yes'; 'Yes' 'No'; 'No' 'No'};
            
            for ii = 1:3
                spikeFile = sprintf('%s\\%s_spikes_concat_forceclustered_%s_ignorenoise_%s.mat',fileDir,filename,fileOrder{ii,1},fileOrder{ii,2});

                if exist(spikeFile,'file')
                    load(spikeFile,'sortedSpikeTimes','sortedSpikeClusters','sortedSpikeChannels');
                end
            end
        end
        
        assert(exist('sortedSpikeTimes','var') > 0,'No spike times file found');

        widths = valuess{opts.confvar}; %#ok<USENS> %*25/2; %#ok<USENS>
        shifts = valuess{opts.infovar}; %*8; TODO : specify this as a function in varargin

        nWidths = levels(opts.confvar);
        nPhases = levels(opts.infovar);

    %     responsiveChannels = [13 2; 17 1; 35 2; 36 1; 37 1; 37 2; 48 2; 58 1; 62 1; 71 2; 73 1];

        if strcmp(opts.significance,'rf')
            recordings = initRecordings;
            if isfield(recording,'rfFileIndex')
                rfRecording = recordings(recording.rfFileIndex);

                significantResponses = false(nCells,1);

                for ii = 1:nCells
                    rfDir = getAnalysisOutputDir(rfRecording);
                    
                    if iscell(cells)
                        channel = cells{cellIndex,1};
                        cluster = cells{cellIndex,2};
                    else
                        channel = char(cells(ii,1:2));
                        cluster = num2str(cells(ii,3));
                    end
                    
                    rfFile = [rfDir '\srrf_' rfRecording.dataFile '_channel_' channel '_cluster_' cluster '_clustered_phototrigger_all_spikes.mat'];

                    if exist(rfFile,'file')
                        load(rfFile);
                        significantResponses(ii) = logical(significant);
                    end
                end
            else
                significantResponses = true(nCells,1);
            end

            responsiveCells = find(significantResponses);
        elseif strcmp(opts.significance,'visual')
            lightResponses = importdata([fileDir '\lightresponses.txt'],'\t');
            lightResponses = [double(num2str(lightResponses(:,1))) lightResponses(:,2)+1];
            responsiveCells = cellIndices(sub2ind(size(cellIndices),lightResponses(:,1),lightResponses(:,2),lightResponses(:,3)));
        elseif strcmp(opts.significance,'debug')
            lightResponses = [double(num2str(51)) 2];
            responsiveCells = cellIndices(sub2ind(size(cellIndices),lightResponses(:,1),lightResponses(:,2),lightResponses(:,3)));
        elseif strcmp(opts.significance,'file') && ischar(opts.responsefile) && exist(opts.responsefile,'file')
            lightResponses = load(opts.responsefile,'responsiveCells');
            lightResponses = lightResponses.responsiveCells;
            lightResponses = [double(num2str(lightResponses(:,1))) lightResponses(:,2)];
%             responsiveCells = cellIndices(sub2ind(size(cellIndices),lightResponses(:,1),lightResponses(:,2),lightResponses(:,3)));
            responsiveCells = [];
            
            for ii = 1:size(lightResponses,1)
                index = find(ismember(cells,lightResponses(ii,:),'rows'));
                
                if ~isempty(index)
                    responsiveCells(end+1) = index; %#ok<AGROW>
                end
            end
        else
            responsiveCells = 1:nCells;
        end

        nResponses = numel(responsiveCells);

        assert(nResponses <= size(rasters,1));

    %     alpha = 0.05/3200;

    %     for ii = 1:nCells
    %         FRs = [];
    %         validBins = cell(nWidths,nPhases);
    %         binWidths = zeros(nWidths,nPhases);

    %         for jj = 1:nWidths
    %             for kk = 1:nPhases
    %                 edges = edgess{ii,jj,kk,1}; %#ok<USENS>
    %                 binWidths(jj,kk) = mean(diff(edges));
    %                 validBins{jj,kk} = edges >= 0 & edges < stimulusLengths(jj,kk,1);
    %                 FRs = [FRs; histograms{ii,jj,kk,1}(validBins{jj,kk})/binWidths(jj,kk)]; %#ok<AGROW,USENS>
    %             end
    %         end

    %         meanFR = mean(FRs);
    %         stdFR = std(FRs);

    %         for jj = 1:nWidths
    %             for kk = 1:nPhases
    %                 maxFR = max(histograms{ii,jj,kk,1}(validBins{jj,kk}))/binWidths(jj,kk);
    %                 minFR = min(histograms{ii,jj,kk,1}(validBins{jj,kk}))/binWidths(jj,kk);
    %                 significant = ztest(maxFR,meanFR,stdFR,alpha) || ztest(minFR,meanFR,stdFR,alpha);
    %                 significantResponses(ii) = significantResponses(ii) + significant;
    %             end
    %         end
    %     end

    %     responsiveCells = find(significantResponses > nWidths*nPhases*0.05);
    %     nResponses = numel(responsiveCells);

        datas = cell(nResponses,nWidths,nPhases,3,2);
        means = zeros(nResponses,nWidths,nPhases,2,2);
        stds = zeros(nResponses,nWidths,nPhases,2,2);
        ns = zeros(nResponses,nWidths,nPhases,2,2);
        maxCount = -Inf;
        maxLatency = -Inf;
        
        infoVars = [];
        confVars = [];
        respVars = cell(1,nResponses);

        %% DATA MUNGING

    %     f = [figure figure figure];
    %     figure;
        alpha = 0.05/(nWidths*nPhases*nCells);
    %     kernel = normpdf(-25:25,0,5);
        t = (-100:100);
        kernel = zeros(201,3);
        kernel(51:151,1) = 1;
        kernel(101:201,2) = 1;
        kernel(:,3) = normpdf(t,0,5);
        kernel(101:201,4) = normpdf(t(101:201)-50,0,5);
        kernel(:,5) = heaviside(t).*(1-exp(-t)).*exp(-t/20);
        kernel = kernel./repmat(sum(kernel),201,1);
        significantResponses = false(nCells,1);
        for ii = 1:nResponses
            cellIndex = responsiveCells(ii);

            if iscell(cells)
                channel = cells{cellIndex,1};
                cluster = cells{cellIndex,2};
            else
                channel = cells(cellIndex,1:2);
                cluster = cells(cellIndex,3);
            end

            if iscell(sortedSpikeTimes)
                spikeTimes = sortedSpikeTimes{cellIndex};
            else
                spikeTimes = sortedSpikeTimes(sortedSpikeChannels(:,1) == channel(1) & sortedSpikeChannels(:,2) == channel(2) & sortedSpikeClusters == cluster); %#ok<NODEF>
            end
            
            spikeTimes = ceil(1000*spikeTimes);
    %         phat = gamfit(1./diff(spikeTimes));
    %         threshlow = gaminv(0.025,phat(1),phat(2));
    %         threshhigh = gaminv(0.975,phat(1),phat(2));
            ifr = zeros(max(spikeTimes),5);
            ifr(spikeTimes,:) = 1;
            
            for jj = 1:5
                ifr(:,jj) = conv(ifr(:,jj),kernel(:,jj),'same');
            end
            
            mu = mean(ifr);
            sigma = std(ifr);
            threshlow = mu-3*sigma;
            threshhigh = mu+3*sigma;
    %         f = [figure figure];
    %         ax = zeros(nWidths,nPhases,2);

    %         allValidLatencies = {[] []};
            for jj = 1:nWidths
                for kk = 1:nPhases

                    raster = rasters{cellIndex,jj,kk,1};
                    stimulusTiming = stimulusTimings{jj,kk,1};
                    nTrials = numel(raster);
                    counts = zeros(nTrials,1,2);
                    latencies = zeros(nTrials,2);
                    % TODO : other measures of latency
                    firstSpikes = zeros(nTrials,1,2);
                    peakTimes = zeros(nTrials,5,2);
                    
    %                 ifrs = {zeros(ceil(1000*max(diff(stimulusTiming(:,1,:)))),numel(raster)) zeros(ceil(1000*max(diff(stimulusTiming(:,2,:)))),numel(raster))};
    %                 for ll = 1:2
    %                     figure(f(ll));
    %                     ax(jj,kk,ll) = subplot(nWidths,nPhases,nPhases*(jj-1)+kk);
    %                     xlim([0 numel(raster)+1]);
    %                     ylim([0 max(stimulusTiming(ll,2,:) - stimulusTiming(ll,1,:))])
    %                     hold on;
    %                 end

                    if ii == 1
                        infoVars = [infoVars; repmat(shifts(kk),nTrials,1)]; %#ok<AGROW>
                        confVars = [confVars; repmat(widths(jj),nTrials,1)]; %#ok<AGROW>
                    end
                    
                    for ll = 1:nTrials
                        stimulusStart = stimulusTiming(1,1,ll);
                        row = raster{ll};

                        for mm = 1:size(stimulusTiming,2)
    %                         axes(ax(jj,kk,mm));

                            validSpikes = find(row > stimulusTiming(1,mm,ll)-stimulusStart & row <= stimulusTiming(2,mm,ll)-stimulusStart);
                            response = row(validSpikes) - (stimulusTiming(1,mm,ll)-stimulusStart);

                            if isempty(response)
                                counts(ll,1,mm) = 0;
                                firstSpikes(ll,1,mm) = Inf;
                                peakTimes(ll,:,mm) = Inf;
                                continue;
                            end

                            counts(ll,1,mm) = numel(response);
    %                         plot(ll*ones(size(response)),response,'LineStyle','none','Marker','.');

%                             if strcmp(opts.latencytype,'firstspike')
                                firstSpikes(ll,1,mm) = 1000*response(1);
% 
%                                 if ~strcmp(opts.significance,'raster')
%                                     continue;
%                                 end
%                             end

                            if strcmp(opts.latencytype,'burst')
                                if numel(response) == 1
                                    latencies(ll,mm) = 1000*response(1);
                                else
                                    [~,~,burstStarts] = burstold(response);

                                    if isempty(burstStarts)
                                        latencies(ll,mm) = 1000*response(1);
                                    else
                                        latencies(ll,mm) = 1000*response(burstStarts(1));
                                    end
                                end
                            end


                            indices = ceil(1000*(row - diff(stimulusTiming([1 mm],1,ll))));
    %                         ifrs = zeros(ceil(1000*max(diff(stimulusTiming(:,mm,:)))),1);
    %                         
    %                         firstSpike = validSpikes(1);
    %                         if firstSpike > 1
    % %                             ifrs{mm}(1:indices(firstSpike)-1,ll) = 1/diff(row(firstSpike-[1 0]));
    %                             ifrs(1:indices(firstSpike)-1) = 1/diff(row(firstSpike-[1 0]));
    %                         end
    %                         
    %                         for nn = 1:numel(validSpikes)-1
    %                             thisSpike = validSpikes(nn);
    %                             nextSpike = validSpikes(nn+1);
    % %                             ifrs{mm}(indices(thisSpike):indices(nextSpike)-1,ll) = 1/diff(row([thisSpike nextSpike]));
    %                             ifrs(indices(thisSpike):indices(nextSpike)-1) = 1/diff(row([thisSpike nextSpike]));
    %                         end
    %                         
    %                         lastSpike = validSpikes(end);
    %                         if lastSpike < numel(row);
    % %                             ifrs{mm}(indices(lastSpike):end,ll) = 1/diff(row(lastSpike+[0 1]));
    %                             ifrs(indices(lastSpike):end) = 1/diff(row(lastSpike+[0 1]));
    %                         end
    %                         
    %                         plow = gamcdf(ifrs,phat(1),phat(2));
    %                         phigh = 1-plow;
    %                         significant = (plow < alpha/2 | phigh < alpha/2);

    %                         indices = ceil(row*1000);
%                             ifrs = zeros(ceil(1000*max(diff(stimulusTiming(:,mm,:))))+2000,1);
%                             indices = indices(indices > -1000 & indices + 1000 < numel(ifrs));
%                             ifrs(indices+1000) = 1;
%                             ifrs = conv(ifrs,kernel,'same');
%                             ifrs = ifrs(1001:end-1000);

                            ifrs = zeros(ceil(1000*max(diff(stimulusTiming(:,mm,:)))),5);
                            indices = indices(indices > 0 & indices <= numel(ifrs));
                            ifrs(indices,:) = 1;
                            
                            for nn = 1:5
                                ifrs(:,nn) = conv(ifrs(:,nn),kernel(:,nn),'same');
                            end

                            if strcmp(opts.significance,'raster')
                                Z = (ifrs-mu)/sigma;
                                p = 2*normcdf(-abs(Z));

                                if any(p < alpha);
                                    significantResponses(ii) = true;
                                end
                            end

    % %                         clf;
    % %                         hold on;
    % %                         plot(ifrs);
    %                         line([0 0; numel(ifrs) numel(ifrs)],[threshlow threshhigh; threshlow threshhigh],'LineStyle','--','Color','k')
    % %                         xlim([0 numel(ifrs)]);#
    %                         
    %                         hold on;
    %                         plot(ifrs);
    % %                         line(repmat([validIndices(1);validIndices(end)],1,2),repmat(mu+2*sigma*[-1 1],2,1),'LineStyle','--','Color','k')
    % %                         xlim([0 numel(ifrs)]);
    %                         fprintf('%d\t%d\t%d\t%d\t%d\n',ii,jj,kk,ll,mm)
    % %                         clf;

                            if strcmp(opts.latencytype,'response')
                                if ifrs(1) < threshlow || ifrs(1) > threshhigh
                                    prevResponseEnd = find(ifrs >= threshlow & ifrs <= threshhigh,1);

                                    latency = find(ifrs(prevResponseEnd+1:end) < threshlow | ifrs(prevResponseEnd+1:end) > threshhigh,1)+prevResponseEnd;
                                else
                                    latency = find(ifrs < threshlow | ifrs > threshhigh,1);
                                end
                                
                                if isempty(latency)
                                    latencies(ll,mm) = Inf;
                                else
                                    latencies(ll,mm) = latency;
                                end
                            end
                            
%                             elseif strcmp(opts.latencytype,'peak')
                                latency = zeros(1,5);
                                
                                for nn = 1:5
                                    latency(nn) = find(ifrs(:,nn) == max(ifrs(:,nn)),1);
                                end
%                             else
%                                 continue;
%                             end

                            if isempty(latency) % can this even happen???
                                peakTimes(ll,:,mm) = Inf;
                            else
    %                             line([ll ll; ll ll] + [-1 -1; 1 1]/2,latency*ones(2,2)/1000,'Color','r');
                                peakTimes(ll,:,mm) = latency;
                            end
    %                         clf;
                    
                        end
                    end
                    
                    resp = [counts firstSpikes peakTimes];
                    counts = squeeze(counts);
                    
                    if isempty(respVars{ii})
                        respVars{ii} = resp;
                    else
                        respVars{ii} = [respVars{ii}; resp];
                    end
                    
                    if strcmp(opts.latencytype,'firstspike')
                        latencies = squeeze(firstSpikes);
                    elseif strcmp(opts.latencytype,'peak')
                        latencies = squeeze(peakTimes(:,:,1));
                    end

                    maxCount = max(maxCount,max(max(counts)));
                    means(ii,jj,kk,1,:) = mean(counts);
                    stds(ii,jj,kk,1,:) = std(counts);
                    ns(ii,jj,kk,1,:) = size(counts,1);

    %                 clf;
                    for ll = 1:2
                        datas{ii,jj,kk,1,ll} = counts(:,ll);
                        validLatencies = latencies(isfinite(latencies(:,ll)),ll);

                        if ~isempty(validLatencies)
                            maxLatency = max(maxLatency,max(validLatencies));
                        end

    %                     figure(f(ll));
    %                     subplot(nWidths,nPhases,nPhases*(jj-1)+kk);
    %                     hist(validLatencies);
    %                     allValidLatencies{ll} = [allValidLatencies{ll}; validLatencies];

                        datas{ii,jj,kk,2,ll} = latencies(:,ll);
                        means(ii,jj,kk,2,ll) = mean(validLatencies);
                        stds(ii,jj,kk,2,ll) = std(validLatencies);
                        ns(ii,jj,kk,2,ll) = numel(validLatencies);

    %                     axes(ax(jj,kk,ll));
    %                     line([0 51],[1 1]*mean(validLatencies)/1000,'Color','g');

    %                     datas{ii,jj,kk,3,ll} = ifrs{ll};

    %                     subplot(2,1,ll);
    %                     hold on;
    %                     plot((1:size(ifrs{ll},1))/1000,pseudolog(ifrs{ll}))
    %                     plot((1:size(ifrs{ll},1))/1000,mean(pseudolog(ifrs{ll}),2),'Color','k','LineWidth',3)
    %                     xlim([0 size(ifrs{ll},1)/1000])
                    end
                end
            end

    %         figure(f(3));
    %         subplot(1,2,1);
    %         hist(allValidLatencies{1},50);
    %         
    %         subplot(1,2,2);
    %         hist(allValidLatencies{2},50);
        end
        
        save(sprintf('%s\\%s_info_data',fileDir,filename),'infoVars','confVars','respVars','-v7.3');

        [widths,widthIndices] = sort(widths);
        [shifts,shiftIndices] = sort(shifts);

        datas = datas(:,widthIndices,shiftIndices,:,:);
        means = means(:,widthIndices,shiftIndices,:,:);
        stds = stds(:,widthIndices,shiftIndices,:,:);
        ns = ns(:,widthIndices,shiftIndices,:,:);
        sortedRasters = rasters(:,widthIndices,shiftIndices,:);
        sortedTimings = stimulusTimings(widthIndices,shiftIndices,:);

        if strcmp(opts.significance,'raster')
            nResponses = sum(significantResponses);

            datas = datas(significantResponses,:,:,:,:);
            means = means(significantResponses,:,:,:,:);
            stds = stds(significantResponses,:,:,:,:);
            ns = ns(significantResponses,:,:,:,:);
        end

        %% BASIC PLOTS

    %     colours = 'rgbm';
    %     
    %     legends = cell(numel(widths),1);
    %     
    %     for ii = 1:numel(widths)
    %         legends{ii} = [num2str(widths(ii)) '{\mu}m Bar'];
    %     end
    %     
    %     variables = {'count' [opts.latencytype 'latency']};
    %     responseTypes = {'on' 'off'};
    %     titleVariables = {'Spike Count' 'First Spike Latency'};
    %     titleResponseTypes = {'Grating' 'Mask'};
    %     ylabels = {'Spike Count','First Spike Latency (ms)'};
    %     
    %     [nRows,nCols] = subplots(nResponses+1);
    %    
    %     for gg = 1:2
    %         for hh = 1:2
    %             figure;
    %             set(gcf,'Position',[0 900 1280 1024]);
    % 
    %             for ii = 1:nResponses
    %                 cellIndex = responsiveCells(ii);
    %                 
    %                 subplot(nRows,nCols,ii);
    %                 hold on;
    % 
    %                 for jj = 1:nWidths
    %                     errorbar(shifts,squeeze(means(ii,jj,:,gg,hh)),2*squeeze(stds(ii,jj,:,gg,hh)./sqrt(ns(ii,jj,:,gg,hh))),'Color',colours(jj));
    %                 end
    % 
    %                 xlim([0 7]);
    % 
    %                 if ii == 9
    %                     xlabel('Shift #');
    %                     ylabel(ylabels{gg});
    %                 end
    % 
    %                 title(['Channel ' char(cells(cellIndex,1:2)) ' cluster ' num2str(cells(cellIndex,3))]);
    %             end
    %     
    %             leg = legend(legends);
    % 
    %             fakePlot = subplot(nRows,nCols,nRows*nCols);
    %             axis off;
    % 
    %             set(leg,'Position',get(fakePlot,'Position'));
    %             
    %             set(gcf, 'PaperPositionMode', 'auto');
    %            
    %             %%
    %             mtit([titleVariables{gg} ' vs bar width and bar position during presentation of ' titleResponseTypes{hh}],'xoff',0,'yoff',0.05);
    %             %%
    %             
    %             saveas(gcf,[fileDir '\' filename '_' variables{gg} '_vs_bar_shift_' responseTypes{hh} '.png']);
    %         end
    %     end

        %% PER-STIMULUS RASTERS
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             fig = figure;
    %             set(fig,'Position',[-1920 -180 1920 1080],'PaperPositionMode','auto');
    %             hold on;
    %             cellIndex = responsiveCells(ii);
    % 
    %             for jj = 1:nWidths
    % 
    %     %             if jj < nResponses
    %     %                 line(600*jj*ones(1,2),[0 400],'LineStyle','--');
    %     %             end
    % 
    %                 for kk = 1:nPhases
    %     %                 if kk < nPhases
    %     %                     line([0 600*nResponses],50*kk*ones(1,2),'LineStyle','--');
    %     %                 end 
    % 
    %                     raster = sortedRasters{cellIndex,jj,kk,1};
    %     %                 subplot(nPhases,nResponses,nResponses*(kk-1)+jj);
    %     %                 subplot('Position',[(jj-1)/nResponses (kk-1)/nPhases 1920/nResponses-10 1080/nPhases-10]);
    % %                     subplot(nWidths,nPhases,nPhases*(jj-1)+kk)
    %                     subplot('Position',[(230*(kk-1)+50)/1920 (250*(4-jj)+50)/1080 210/1920 230/1080])
    %                     set(gca,'YDir','reverse');
    %                     
    %                     if jj == 1
    %                         title(sprintf('Phase %4.3f',(kk-1)/nPhases));
    %                     end
    %                     
    %                     if jj < 4
    %                         set(gca,'XTick',[]);
    %                     elseif kk == 1
    %                         xlabel('Trial #');
    %                         ylabel('Time from grating onset (s)');
    %                     end
    %                     
    %                     if kk > 1
    %                         set(gca,'YTick',[]);
    %                         
    %                         if kk == nPhases
    %                             set(gca,'YAxisLocation','right');
    %                             yLabel = ylabel(['Width ' num2str(widths(jj)) '{\mu}m']);
    %                             set(yLabel,'Rotation',270);
    %                             set(yLabel,'VerticalAlignment','bottom')
    %                         end
    %                     end
    %                     
    %                     hold on;
    % 
    %                     stimulusTiming = sortedTimings{jj,kk,1};
    %                     xlim([0 numel(raster)]);
    %                     ylim([0 median(stimulusTiming(2,2,:)-stimulusTiming(1,1,:))]);
    % 
    %                     for ll = 1:numel(raster)
    %                         for mm = 1:2
    %                             colour = 0.9*ones(1,3);
    %                             colour(mm) = 1.0;
    %                             fill(ll+[-1 -1 1 1],stimulusTiming([1 2 2 1],mm,ll)-stimulusTiming(1,1,ll),colour,'EdgeColor','none');
    %                         end
    %                     end
    %                         
    %                     for ll = 1:numel(raster)
    %                         row = raster{ll};
    % 
    %                         row = row(row > 0 & row <= stimulusTiming(2,2,ll) - stimulusTiming(1,1,ll));
    %                         
    %                         if isempty(row)
    %                             continue;
    %                         elseif true || numel(row) == 1
    %                             plot(ll*ones(size(row))-0.5,row,'Color','b','LineStyle','none','Marker','.');
    %                             continue;
    %                         end
    %                         
    %                         [~,burstLen,burstStart] = burstold(row);
    %                         
    %                         if isempty(burstLen)
    %                             burstLen = 0;
    %                             burstStart = 1;
    %                         end
    %                         
    %                         n = numel(row);
    %                         isInBurst = (1:n >= burstStart(1) & 1:n < burstStart(1) + burstLen(1));
    %                         
    %                         burst = row(isInBurst);
    %                         notBurst = row(~isInBurst);
    %                         
    %                         plot(ll*ones(size(burst))-0.5,burst,'Color','r','LineStyle','none','Marker','.');
    %                         plot(ll*ones(size(notBurst))-0.5,notBurst,'Color','b','LineStyle','none','Marker','.');
    %                     end
    %                 end
    %             end
    %             
    %             saveas(fig,[fileDir '\' filename '_raster_channel_' char(cells(cellIndex,1:2)) '_cluster_' num2str(cells(cellIndex,3)) '.png']);
    %         end
    %     end

    %     return;

        %% MUTUAL INFORMATION BETWEEN COUNT AND STIMULUS - new method
    %     totalCountEntropy = zeros(nResponses,2);
    %     bins = 0:maxCount;
    %     countPDF = zeros(numel(bins),nResponses,2);
    %     countCPDFs = zeros(numel(bins),nResponses,nWidths,nPhases,2);
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             cellCounts = squeeze(datas(ii,:,:,1,hh));
    %             allCounts = [];
    %             
    %             for jj = 1:nWidths
    %                 for kk = 1:nPhases
    %                     stimulusCounts = cellCounts{jj,kk};
    %                     countCPDFs(:,ii,jj,kk,hh) = hist(stimulusCounts,bins)/numel(stimulusCounts);
    %                     allCounts = [allCounts; stimulusCounts]; %#ok<AGROW>
    %                 end
    %             end
    %             
    %             pdf = hist(allCounts,bins)/numel(allCounts);
    %             countPDF(:,ii,hh) = pdf;
    %             informationContext = pdf.*log2(pdf);
    %             informationContext(isnan(informationContext)) = 0; % 0*log(0) = 0
    %             totalCountEntropy(ii,hh) = -sum(informationContext);
    %             
    %             % using formula H(Y|X) = sum_x(p(x)*sum_y(p(y|x)*log(1/p(y|x)))),
    %             % noting log(1/x) = -log(x) and all values of X (here stimulus)
    %             % are equiprobable, hence sum_x(p(x)*f(x)) = mean(f(x))
    % %             conditionalInformationContext = cpdfs.*(-log2(cpdfs));
    % %             conditionalInformationContext(isnan(conditionalInformationContext)) = 0;
    % %             conditionalCountEntropy(ii,hh) = mean(sum(conditionalInformationContext));
    %         end
    %     end
    %     
    %     countWidthCPDFs = mean(countCPDFs,4);
    %     countWidthICs = countWidthCPDFs.*(-log2(countWidthCPDFs));
    %     countWidthICs(isnan(countWidthICs)) = 0;
    %     countWidthEntropy = squeeze(mean(sum(countWidthICs,1),3));
    %         
    %     countPhaseCPDFs = mean(countCPDFs,3);
    %     countPhaseICs = countPhaseCPDFs.*(-log2(countPhaseCPDFs));
    %     countPhaseICs(isnan(countPhaseICs)) = 0;
    %     countPhaseEntropy = squeeze(mean(sum(countPhaseICs,1),4));
    %     
    %     countStimulusICs = countCPDFs.*(-log2(countCPDFs));
    %     countStimulusICs(isnan(countStimulusICs)) = 0;
    %     countStimulusEntropy = squeeze(mean(mean(sum(countStimulusICs,1),3),4));
    %     
    %     mutualCountWidthInformation = totalCountEntropy - countWidthEntropy;
    %     mutualCountPhaseInformation = totalCountEntropy - countPhaseEntropy;
    %     mutualCountStimulusInformation = totalCountEntropy - countStimulusEntropy;

        %% MUTUAL COUNT INFORMATION DIRECTLY FROM JOINT PDFs
        %* gives the same results as calculating based on the mutual entropy to
        %* within an error of +/- 1e-14

    %     bins = 0:maxCount;
    %     JPDF = zeros(numel(bins),nWidths,nPhases,nResponses,2);
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             for jj = 1:4
    %                 for kk = 1:nPhases
    %                     counts = datas{ii,jj,kk,1,hh};
    %                     JPDF(:,jj,kk,ii,hh) = hist(counts,bins)/numel(counts);
    %                 end
    %             end
    %         end
    %     end
    %     
    %     JPDF = JPDF/(nWidths*nPhases);
    %     
    %     countPDF = sum(sum(JPDF,2),3);
    %     
    %     mutualCountPhaseInformationPerWidth = zeros(4,nResponses,2);
    %     
    %     for ii = 1:nWidths
    %         countPhaseJPDF = nWidths*JPDF(:,ii,:,:,:);
    %         phasePDF = sum(countPhaseJPDF);
    %         countPhasePDFProduct = zeros(size(countPhaseJPDF));
    %         
    %         for jj = 1:2
    %             for kk = 1:nResponses
    %                 countPhasePDFProduct(:,1,:,kk,jj) = squeeze(countPDF(:,1,1,kk,jj))*squeeze(phasePDF(1,1,:,kk,jj))';
    %             end
    %         end
    %         
    %         mutualCountPhaseInformationPerWidth(ii,:,:) = squeeze(sum(sum(plogp(countPhaseJPDF,countPhaseJPDF./countPhasePDFProduct))));
    %     end
    %     
    %     totalCountEntropy2 = -sum(plogp(countPDF));
    %     
    %     countWidthJPDF = sum(JPDF,3);
    %     countPhaseJPDF = sum(JPDF,2);
    %     countStimulusJPDF = reshape(JPDF,[numel(bins) 32 nResponses 2]);
    %     
    %     widthPDF = sum(countWidthJPDF);
    %     phasePDF = sum(countPhaseJPDF);
    %     stimulusPDF = sum(countStimulusJPDF);
    %     
    %     countWidthPDFProduct = zeros(size(countWidthJPDF));
    %     countPhasePDFProduct = zeros(size(countPhaseJPDF));
    %     countStimulusPDFProduct = zeros(size(countStimulusJPDF));
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             countWidthPDFProduct(:,:,1,ii,hh) = squeeze(countPDF(:,1,1,ii,hh))*squeeze(widthPDF(1,:,1,ii,hh));
    %             countPhasePDFProduct(:,1,:,ii,hh) = squeeze(countPDF(:,1,1,ii,hh))*squeeze(phasePDF(1,1,:,ii,hh))';
    %             countStimulusPDFProduct(:,:,ii,hh) = squeeze(countPDF(:,1,1,ii,hh))*squeeze(stimulusPDF(1,:,ii,hh));
    %         end
    %     end
    %     
    %     mutualCountWidthInformation = squeeze(sum(sum(plogp(countWidthJPDF,countWidthJPDF./countWidthPDFProduct))));
    %     mutualCountPhaseInformation = squeeze(sum(sum(plogp(countPhaseJPDF,countPhaseJPDF./countPhasePDFProduct))));
    %     mutualCountStimulusInformation = squeeze(sum(sum(plogp(countStimulusJPDF,countStimulusJPDF./countStimulusPDFProduct))));

        %% MUTUAL INFORMATION BETWEEN COUNT AND STIMULUS AS A FUNCTION OF PHASE

        maxSplit = 6;
        repeats = 10;
        nSplits = 1+repeats*sum(2:maxSplit);
        bins = 0:maxCount;

        splits = [1;kron([2;2;3;3;3;4;4;4;4;5;5;5;5;5;6;6;6;6;6;6],ones(repeats,1))];
        unbiasedCountInfo = zeros(nResponses,nWidths,2);
        
        for ii = 1:nResponses
            tic;
            countPDF = zeros(numel(bins),nWidths,2,2,nSplits);
            phasePDF = zeros(nPhases,nWidths,2,2,nSplits);
            jointPDF = zeros(numel(bins),nPhases,nWidths,2,2,nSplits);
            mutualCountPhaseInformationPerWidth2 = zeros(nWidths,2,2,nSplits);

            for hh = 1:2
                for jj = 1:nWidths
                    allCounts = cell(nSplits,1);

                    for kk = 1:nPhases
                        counts = datas{ii,jj,kk,1,hh};
    %                     allCounts = [allCounts; counts]; %#ok<AGROW>
                        n = numel(counts);

                        nn = 0;
                        for ll = 1:maxSplit
                            m = floor(n/ll);

                            for oo = 1:(ll == 1) + (ll > 1)*repeats
                                index = randperm(n);

                                for mm = 1:ll
                                    nn = nn + 1;

                                    if mm == ll
                                        subCounts = counts(index((mm-1)*m+1:n));
                                    else
                                        subCounts = counts(index((mm-1)*m+1:mm*m));
                                    end

                                    allCounts{nn} = [allCounts{nn}; subCounts];

                                    jointPDF(:,kk,jj,hh,1,nn) = hist(subCounts,bins)/numel(subCounts);

                                    if ~any(counts > 0)
                                        jointPDF(:,kk,jj,hh,2,nn) = 0;
                                    else
                                        jointPDF(:,kk,jj,hh,2,nn) = [nan; hist(subCounts(subCounts > 0),bins(2:end))'/numel(subCounts(subCounts > 0))];
                                    end

                                    phasePDF(kk,jj,hh,:,nn) = numel(subCounts);
                                end
                            end
                        end
                    end

                    for kk = 1:nSplits
                        countPDF(:,jj,hh,1,kk) = hist(allCounts{kk},bins)/numel(allCounts{kk});

                        if ~any(allCounts{kk} > 0)
                            countPDF(:,jj,hh,2,kk) = 0;
                        else
                            allCount = allCounts{kk};
                            countPDF(:,jj,hh,2,kk) = [nan; hist(allCount(allCount > 0),bins(2:end))'/numel(allCount(allCount > 0))];
                        end

                        for gg = 1:2
                            phasePDF(:,jj,hh,gg,kk) = phasePDF(:,jj,hh,gg,kk)/sum(phasePDF(:,jj,hh,gg,kk));
                        end
                    end
                end
            end

            jointPDF = jointPDF/nPhases;

            for ff = 1:nSplits
                for hh = 1:2
                    for jj = 1:nWidths
                        for gg = 1:2
                            px = squeeze(countPDF(gg:end,jj,hh,gg,ff));
                            py = squeeze(phasePDF(:,jj,hh,gg,ff));
                            pxpy = repmat(px,1,nPhases).*repmat(py',numel(bins)-gg+1,1);
                            pxy = squeeze(jointPDF(gg:end,:,jj,hh,gg,ff));
                            mutualCountPhaseInformationPerWidth2(jj,hh,gg,ff) = sum(sum(plogp(pxy,pxy./pxpy)));
                        end
                    end
                end
            end

    %     figure;
    %         clf;

            for jj = 1:nWidths
                for kk = 1:2
                    info1 = squeeze(mutualCountPhaseInformationPerWidth2(jj,kk,2,:));
    %                 subplot(2,4,4*(kk-1)+jj);
    %                 hold on;
    %                 plot(splits,info1,'Color','b','LineStyle','none','Marker','o');
                    info2 = zeros(maxSplit,1);
                    info2(1) = info1(1);

                    ll = 1;
                    for mm = 2:maxSplit
                        info2(mm) = mean(info1(ll+(1:repeats*mm)));
                        ll = ll + repeats*mm;
                    end

    %                 plot(1:maxSplit,info2,'Color','g');

                    p = polyfit(splits,info1,2);

    %                 fplot(@(x) polyval(p,x),[0 maxSplit],'Color','r');

                    unbiasedCountInfo(ii,jj,kk) = polyval(p,0);
                end
            end
            toc;
        end
        
        if ee == 1
            unbiasedCountInfos = zeros([size(unbiasedCountInfo) nRecordings]);
            unbiasedCountInfos(:,:,:,1) = unbiasedCountInfo;
            indices = 1:size(unbiasedCountInfos,1);
        else
            indices = [];
            for ii = 1:nResponses
                if isstruct(cellIndices)
                    indices(end+1) = cellIndices.(cells{ii,1})(cells{ii,2}); %#ok<AGROW>
                else
                    indices(end+1) = cellIndices(cells(ii,1),cells(ii,2),cells(ii,3)+1); %#ok<AGROW>
                end
            end
            
            if max(indices) > size(unbiasedCountInfos,1)
                unbiasedCountInfos(end+1:max(indices),:,:,:) = 0;
            end
                
            unbiasedCountInfos(indices,:,:,ee) = unbiasedCountInfo;
        end

        %% MUTUAL INFORMATION BETWEEN COUNT AND SPATIAL FREQUENCY

    %     bins = 0:maxCount;
    %     countPDF = zeros(numel(bins),nResponses,2,2);
    %     widthPDF = zeros(nWidths,nResponses,2,2);
    %     jointPDF = zeros(numel(bins),nWidths,nResponses,2,2);
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             allCounts = [];
    %             
    %             for jj = 1:nWidths
    %                 widthCounts = [];
    %             
    %                 for kk = 1:nPhases
    %                     counts = datas{ii,jj,kk,1,hh};
    %                     widthCounts = [widthCounts; counts]; %#ok<AGROW>
    %                 end
    %                 
    %                 jointPDF(:,jj,ii,hh,1) = hist(widthCounts,bins)/numel(widthCounts);
    % 
    %                 if ~any(widthCounts > 0)
    %                     jointPDF(:,jj,ii,hh,2) = 0;
    %                 else
    %                     jointPDF(:,jj,ii,hh,2) = [nan; hist(widthCounts(widthCounts > 0),bins(2:end))'/numel(widthCounts(widthCounts > 0))];
    %                 end
    %                     
    %                 widthPDF(jj,ii,hh,:) = numel(widthCounts);
    %                 
    %                 allCounts = [allCounts; widthCounts];
    %             end
    %                 
    %             countPDF(:,ii,hh,1) = hist(allCounts,bins)/numel(allCounts);
    %                 
    %             if ~any(allCounts > 0)
    %                 countPDF(:,ii,hh,2) = 0;
    %             else
    %                 countPDF(:,ii,hh,2) = [nan; hist(allCounts(allCounts > 0),bins(2:end))'/numel(allCounts(allCounts > 0))];
    %             end
    %                 
    %             for gg = 1:2
    %                 widthPDF(:,ii,hh,gg) = widthPDF(:,ii,hh,gg)/sum(widthPDF(:,ii,hh,gg));
    %             end
    %         end
    %     end
    %     
    %     jointPDF = jointPDF/4nWidths
    %     mutualCountWidthInformationPerPhase = zeros(nResponses,2,2);
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             for gg = 1:2
    %                 px = squeeze(countPDF(gg:end,ii,hh,gg));
    %                 py = squeeze(widthPDF(:,ii,hh,gg));
    %                 pxpy = repmat(px,1,nWidths).*repmat(py',numel(bins)-gg+1,1);
    %                 pxy = squeeze(jointPDF(gg:end,:,ii,hh,gg));
    %                 mutualCountWidthInformationPerPhase(ii,hh,gg) = sum(sum(plogp(pxy,pxy./pxpy)));
    %             end
    %         end
    %     end

        %% MUTUAL INFORMATION BETWEEN COUNT AND STIMULUS - Moddemeijer method
        %* gives nonsense values but I assume that's to do with using a method
        %* designed for continuous RVs on discrete RVs

    %     totalCountEntropy = zeros(nResponses,2);
    %     conditionalCountWidthEntropy = zeros(nResponses,nWidths,2);
    %     conditionalCountPhaseEntropy = zeros(nResponses,nPhases,2);
    %     conditionalCountStimulusEntropy = zeros(nResponses,nWidths*nPhases,2);
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             allCounts = [];
    %             widthCounts = cell(nWidths,1);
    %             phaseCounts = cell(nPhases,1);
    %             
    %             for jj = 1:nWidths
    %                 for kk = 1:nPhases
    %                     counts = datas{ii,jj,kk,1,hh};
    %                     conditionalCountStimulusEntropy(ii,sub2ind([nWidths,nPhases],jj,kk),hh) = entropy(counts',[],'unbiased',2);
    %                 	widthCounts{jj} = [widthCounts{jj}; counts];
    %                     phaseCounts{kk} = [phaseCounts{kk}; counts];
    %                     allCounts = [allCounts; counts]; %#ok<AGROW>
    %                 end
    %             end
    %             
    %             for jj = 1:nWidths
    %                 conditionalCountWidthEntropy(ii,jj,hh) = entropy(widthCounts{jj}',[],'unbiased',2);
    %             end
    %             
    %             for jj = 1:nPhases
    %                 conditionalCountPhaseEntropy(ii,jj,hh) = entropy(phaseCounts{jj}',[],'unbiased',2);
    %             end
    %             
    %             totalCountEntropy(ii,hh) = entropy(allCounts',[],'unbiased',2);
    %         end
    %     end
    %     
    %     mutualCountWidthInformation3 = totalCountEntropy - squeeze(mean(conditionalCountWidthEntropy,2));
    %     mutualCountPhaseInformation3 = totalCountEntropy - squeeze(mean(conditionalCountPhaseEntropy,2));
    %     mutualCountStimulusInformation3 = totalCountEntropy - squeeze(mean(conditionalCountStimulusEntropy,2));

        %% MUTUAL INFORMATION BETWEEN COUNT AND STIMULUS - old method
        %* |Iold(X;Y) - Inew(X;Y)| < 1e-13
    %     
    %     totalCountEntropy = zeros(nResponses,2);
    %     conditionalCountEntropy = zeros(nResponses,nWidths,nPhases,2);
    %     bins = 0:maxCount;
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             totalCountProbability = zeros(size(bins));
    %             
    %             for jj = 1:nWidths
    %                 for kk = 1:nPhases
    %                     conditionalCountProbability = hist(datas{ii,jj,kk,1,hh},bins);
    %                     totalCountProbability = totalCountProbability + conditionalCountProbability;
    %                     conditionalCountProbability = conditionalCountProbability/ns(ii,jj,kk,1,hh);
    %                     valid = find(conditionalCountProbability > 0);
    %                     conditionalCountEntropy(ii,jj,kk,hh) = -sum(conditionalCountProbability(valid).*log(conditionalCountProbability(valid))/log(2));
    %                 end
    %             end
    %             
    %             totalCountProbability = totalCountProbability/sum(sum(ns(ii,:,:,1,hh),2),3);
    %             valid = find(totalCountProbability > 0);
    %             totalCountEntropy(ii,hh) = -sum(totalCountProbability(valid).*log(totalCountProbability(valid))/log(2));
    %         end
    %     end
    %     
    %     mutualCountInformation = repmat(reshape(totalCountEntropy,[nResponses 1 1 2]),[1 nWidths nPhases 1]) - conditionalCountEntropy;

        %% MUTUAL INFORMATION BETWEEN LATENCY AND STIMULUS - new method
    %     totalLatencyEntropy = zeros(nResponses,2);
    %     conditionalLatencyEntropy = zeros(nResponses,2);
    %     bins = 0:maxLatency;
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             cellLatencies = squeeze(datas(ii,:,:,2,hh));
    %             cpdfs = zeros(numel(bins),nWidths*nPhases);
    %             chInfs = zeros(1,32);
    %             allLatencies = [];
    %             
    %             for jj = 1:nWidths
    %                 for kk = 1:nPhases
    %                     stimulusLatencies = cellLatencies{jj,kk};
    %                     stimulusIndex = sub2ind([nWidths nPhases],jj,kk);
    %                     validLatencies = stimulusLatencies(isfinite(stimulusLatencies));
    %                     
    %                     if any(isinf(stimulusLatencies))
    %                         cpInf = sum(isinf(stimulusLatencies))/numel(stimulusLatencies);
    %                         chInfs(stimulusIndex) = cpInf*(-log2(cpInf));
    %                     end
    %                     
    %                     cpdfs(:,stimulusIndex) = hist(validLatencies,bins)/numel(stimulusLatencies);
    %                     allLatencies = [allLatencies; stimulusLatencies]; %#ok<AGROW>
    %                 end
    %             end
    %             
    %             validLatencies = allLatencies(isfinite(allLatencies));
    %             pdf = hist(validLatencies,bins)/numel(allLatencies);
    %             informationContext = pdf.*log2(pdf);
    %             informationContext(isnan(informationContext)) = 0;
    %             
    %             if any(isinf(allLatencies))
    %                 pInf = sum(isinf(allLatencies))/numel(allLatencies);
    %                 hInf = -pInf*log2(pInf);
    %             else
    %                 hInf = 0;
    %             end
    %             
    %             totalLatencyEntropy(ii,hh) = -trapz(informationContext) + hInf;
    %             
    %             % Exactly the same logic as for count information except
    %             % latency is a continuous variable so we're integrating over
    %             % the discrete estimate of the information context rather than
    %             % summing.  Stimulus is still a discrete variable so mean is
    %             % still appropriate.
    %             conditionalInformationContext = cpdfs.*(-log2(cpdfs));
    %             conditionalInformationContext(isnan(conditionalInformationContext)) = 0;
    %             conditionalLatencyEntropy(ii,hh) = mean(trapz(conditionalInformationContext) + chInfs);
    %         end
    %     end
    %     
    %     mutualLatencyInformation = totalLatencyEntropy - conditionalLatencyEntropy;

        %% MUTUAL INFORMATION BETWEEN LATENCY AND STIMULUS - Moddemeijer method
        %* Entropy of all latencies is a direct calculation using Moddemeijer's
        %* method.  Conditional entropy is computed using the standard formula 
        %* H(Y|X) = sum_x(p(x)H(Y|X = x)), which reduces to mean_x(H(Y|X = x))
        %* because all X are equiprobable, with the inner part computed by 
        %* applying Moddemeijer's method to the latencies in each condition.
        %* However, we have to make a slight modification to account for
        %* infinite latencies, but following Gollisch & Meister (2008), this
        %* amounts to adding a delta-like function with area equal to the
        %* proportion of non-responses (i.e. infinite latencies).  The finite
        %* part of the entropy integrals is computed using Moddemeijer's
        %* method and infinite part computed separately and added on at the end.
    %     totalLatencyEntropy = zeros(nResponses,2);
    %     conditionalLatencyWidthEntropy = zeros(nResponses,nWidths,2);
    %     conditionalLatencyPhaseEntropy = zeros(nResponses,nPhases,2);
    %     conditionalLatencyStimulusEntropy = zeros(nResponses,nWidths*nPhases,2);
    %     conditionalLatencyPhaseEntropyPerWidth = zeros(nResponses,nPhases,nWidths,2);
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             allValidLatencies = [];
    %             widthValidLatencies = cell(nWidths,1);
    %             phaseValidLatencies = cell(nPhases,1);
    % %             phaseValidLatenciesPerWidth = cell(nPhases,nWidths);
    %             nInfAll = 0;
    %             nInfWidth = zeros(nWidths,1);
    %             nInfPhase = zeros(nPhases,1);
    % %             nInfPhasePerWidth = zeros(nPhases,nWidths);
    %             
    %             for jj = 1:nWidths
    %                 for kk = 1:nPhases
    %                     latencies = datas{ii,jj,kk,2,hh};
    %                     validLatencies = latencies(isfinite(latencies));
    %                     
    %                     if isempty(validLatencies)
    %                         h = 0;
    %                     else
    %                         h = entropy(validLatencies',[],'unbiased',2);
    %                     end
    %                     
    %                 	widthValidLatencies{jj} = [widthValidLatencies{jj}; validLatencies];
    %                     phaseValidLatencies{kk} = [phaseValidLatencies{kk}; validLatencies];
    % %                     phaseValidLatenciesPerWidth{kk,jj} = validLatencies;
    %                     allValidLatencies = [allValidLatencies; validLatencies]; %#ok<AGROW>
    %                     infLatencies = latencies(isinf(latencies));
    %                     
    %                     if ~isempty(infLatencies)
    %                         nInf = numel(infLatencies);
    %                         pInf = nInf/numel(latencies);
    %                         h = h - pInf*log2(pInf);
    %                         nInfWidth(jj) = nInfWidth(jj) + nInf;
    %                         nInfPhase(kk) = nInfPhase(kk) + nInf;
    % %                         nInfPhasePerWidth(kk,jj) = nInf;
    %                         nInfAll = nInfAll + nInf;
    %                     end
    %                     
    %                     conditionalLatencyStimulusEntropy(ii,sub2ind([nWidths,nPhases],jj,kk),hh) = h;
    %                     conditionalLatencyPhaseEntropyPerWidth(ii,kk,jj,hh) = h;
    %                 end
    %             end
    %             
    %             for jj = 1:nWidths
    %                 if isempty(widthValidLatencies{jj})
    %                     h = 0;
    %                 else
    %                     h = entropy(widthValidLatencies{jj}',[],'unbiased',2);
    %                 end
    %                 
    %                 if nInfWidth(jj) > 0
    %                     pInfWidth = nInfWidth(jj)/(nInfWidth(jj) + numel(widthValidLatencies{jj}));
    %                     h = h - pInfWidth*log2(pInfWidth);
    %                 end
    %                 
    %                 conditionalLatencyWidthEntropy(ii,jj,hh) = h;
    %             end
    %             
    %             for jj = 1:nPhases
    %                 if isempty(phaseValidLatencies{jj})
    %                     h = 0;
    %                 else
    %                     h = entropy(phaseValidLatencies{jj}',[],'unbiased',2);
    %                 end
    %                 
    %                 if nInfPhase(jj) > 0
    %                     pInfPhase = nInfPhase(jj)/(nInfPhase(jj) + numel(phaseValidLatencies{jj}));
    %                     h = h - pInfPhase*log2(pInfPhase);
    %                 end
    %                 
    %                 conditionalLatencyPhaseEntropy(ii,jj,hh) = h;
    %             end
    %             
    % %             for jj = 1:nWidths
    % %                 for kk = 1:nPhases
    % %                     if isempty(phaseValidLatenciesPerWidth{kk,jj})
    % %                         h = 0;
    % %                     else
    % %                         h = entropy(phaseValidLatenciesPerWidth{kk,jj});
    % %                     end
    % %                     
    % %                     if nInfPhasePerWidth(kk,jj) > 0
    % %                         pInfPhasePerWidth = nInfPhasePerWidth(kk,jj)/(nInfPhasePerWidth(kk,jj) + numel(phaseValidLatenciesPerWidth{kk,jj}));
    % %                         h = h - pInfPhasePerWidth*log2(pInfPhasePerWidth);
    % %                     end
    % %                     
    % %                     conditionalLatencyPhaseEntropyPerWidth(ii,kk,jj,hh) = h;
    % %                 end
    % %             end
    %             
    %             if isempty(allValidLatencies)
    %                 h = 0;
    %             else
    %                 h = entropy(allValidLatencies',[],'unbiased',2);
    %             end
    %             
    %             if nInfAll > 0
    %                 pInfAll = nInfAll/(nInfAll + numel(allValidLatencies));
    %                 h = h - pInfAll*log2(pInfAll);
    %             end
    %                 
    %             totalLatencyEntropy(ii,hh) = h;
    %         end
    %     end
    %     
    %     mutualLatencyWidthInformation = totalLatencyEntropy - squeeze(mean(conditionalLatencyWidthEntropy,2));
    %     mutualLatencyPhaseInformation = totalLatencyEntropy - squeeze(mean(conditionalLatencyPhaseEntropy,2));
    %     mutualLatencyStimulusInformation = totalLatencyEntropy - squeeze(mean(conditionalLatencyStimulusEntropy,2));
    %     
    %     % conditionalLatencyWidthEntropy, before taking the mean, is just the
    %     % total latency entropy for each cell considering each grating width
    %     % separately
    %     mutualLatencyPhaseInformationPerWidth = conditionalLatencyWidthEntropy - squeeze(mean(conditionalLatencyPhaseEntropyPerWidth,2));

        %% MUTUAL INFORMATION BETWEEN LATENCY AND STIMULUS PER WIDTH - Moddemeijer method
    %     totalLatencyEntropy = zeros(nResponses,nWidths,2);
    %     conditionalLatencyPhaseEntropy = zeros(nResponses,nPhases,nWidths,2);
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             for jj = 1:nWidths
    %                 allLatencies = [];
    %                 
    % %                 close all;
    %                 
    %                 for kk = 1:nPhases
    %                     latencies = datas{ii,jj,kk,2,hh};
    %                     conditionalLatencyPhaseEntropy(ii,kk,jj,hh) = moddemeijerEntropy(latencies);
    %                     allLatencies = [allLatencies; latencies]; %#ok<AGROW>
    %                     
    % %                     valid = latencies(isfinite(latencies));
    % %                     
    % %                     if numel(valid) < 3;
    % %                         continue;
    % %                     end
    % %                     
    % %                     n = ceil(numel(valid)/2)-1;
    % %                     g = zeros(n,1);
    % %                     for ll = 1:n
    % %                         g(ll) = vasicekEntropy(valid,ll,true);
    % %                     end
    % %                     
    % %                     h = entropy(valid');
    % %                     
    % %                     figure;
    % %                     
    % %                     plot(1:n,g,1:n,h*ones(size(g)));
    %                 end
    %                 
    %                 totalLatencyEntropy(ii,jj,hh) = moddemeijerEntropy(allLatencies);
    %                 
    % %                 allValid = allLatencies(isfinite(allLatencies));
    % %                 
    % %                 n = ceil(numel(allValid)-1);
    % %                 g = zeros(n,1);
    % %                 
    % %                 for ll = 1:n
    % %                     g(ll) = vasicekEntropy(allValid,ll,true);
    % %                 end
    % %                     
    % %                 figure;
    % % 
    % %                 plot(1:n,g,1:n,entropy(allValid')*ones(size(g)));
    %             end
    %         end
    %     end
    %     
    %     mutualLatencyPhaseInformationPerWidth2 = totalLatencyEntropy - squeeze(mean(conditionalLatencyPhaseEntropy,2));

        %% MUTUAL INFORMATION BETWEEN LATENCY AND STIMULUS PER WIDTH - Moddemeijer method plus updated correction for delta at infinity

    %     totalLatencyEntropy2 = zeros(nResponses,nWidths,2,4);
    %     conditionalLatencyPhaseEntropy2 = zeros(nResponses,nWidths,2,4);
        mutualLatencyPhaseInformationPerWidth4 = zeros(nResponses,nWidths,2,nSplits);
    %     mutualInverseLatencyPhaseInformationPerWidth4 = zeros(nResponses,nWidths,2);

        for hh = 1:2
            for ii = 1:nResponses
                for jj = 1:nWidths
    %                 allLatencies = [];
    %                 allPhases = [];
                    allValidLatencies = [];
                    allValidPhases = [];
    %                 conditionalEntropy = nan(nPhases,4);

    %                 close all;

                    for kk = 1:nPhases
                        latencies = datas{ii,jj,kk,2,hh};

                        validLatencies = latencies(isfinite(latencies));
    %                     n = max(ceil(numel(validLatencies)/2)-1,1);

    %                     if ~isempty(validLatencies)
    %                         conditionalEntropy(kk,2) = entropy(validLatencies',[],'unbiased',2);
    %                         conditionalEntropy(kk,4) = entropy(1./validLatencies',[],'unbiased',2);
    %                         conditionalEntropy(kk,2) = nan;
    %                     elseif numel(validLatencies) == 1
    %                         conditionalEntropy(kk,2) = 0;
    %                     else
    %                         g = zeros(n,1);
    %                         
    %                         for ll = 1:n
    %                             g(ll) = vasicekEntropy(validLatencies,ll,true,2);
    %                         end
    %                         
    %                         conditionalEntropy(kk,2) = max(g(isfinite(g)));
    %                     end
    %                     
    %                     h = zeros(n,1);
    %                     
    %                     for ll = 1:n
    %                         h(ll) = infCorrectedEntropy(latencies,@(x) vasicekEntropy(x,ll,true,2),2);
    %                     end
    %                     
    %                     figure;
    %                     plot(1:n,h,1:n,infCorrectedEntropy(latencies,@(x) entropy(x',[],'unbiased',2),2)*ones(size(h)));
    %                     
    % %                     fprintf('%d\t%d\t%d\t%d\n',ii,kk,jj,hh);
    % 
    %                     if isempty(h)
    %                         h = 0;
    %                     end
    %                     
    %                     conditionalEntropy(kk,1) = max(h(isfinite(h)));
    %                     end

    %                     conditionalEntropy(kk,1) = infCorrectedEntropy(latencies,@(x) entropy(x',[],'unbiased',2),2);
    %                     conditionalEntropy(kk,3) = infCorrectedEntropy(1./latencies,@(x) entropy(x',[],'unbiased',2),2);

    %                     allLatencies = [allLatencies; latencies]; %#ok<AGROW>
    %                     allPhases = [allPhases; kk*ones(size(latencies))];  %#ok<AGROW>
                        allValidLatencies = [allValidLatencies; validLatencies]; %#ok<AGROW>
                        allValidPhases = [allValidPhases; kk*ones(size(validLatencies))]; %#ok<AGROW>
                    end

    %                 for kk = 1:nWidths
    %                     h = conditionalEntropy(:,kk);
    %                     h = h(isfinite(h));
    %                     
    %                     if isempty(h)
    %                         h = 0;
    %                     else
    %                         h = mean(h);
    %                     end
    %                     
    %                     conditionalLatencyPhaseEntropy2(ii,jj,hh,kk) = h;
    %                 end 

    %                 validLatencies = allLatencies(isfinite(allLatencies));
    %                 n = ceil(numel(validLatencies)/2)-1;

    %                 if ~isempty(validLatencies)
    %                     totalLatencyEntropy2(ii,jj,hh,2) = 0;
    %                     totalLatencyEntropy2(ii,jj,hh,2) = entropy(validLatencies',[],'unbiased',2);
    %                     totalLatencyEntropy2(ii,jj,hh,4) = entropy(1./validLatencies',[],'unbiased',2);
    %                 else
    %                     g = zeros(n,1);
    % 
    %                     for ll = 1:n
    %                         g(ll) = vasicekEntropy(validLatencies,ll,true,2);
    %                     end
    % 
    %                     totalLatencyEntropy2(ii,jj,hh,2) = min(g(isfinite(g)));
    %                 end
    % 
    %                 h = zeros(n,1);
    % 
    %                 for ll = 1:n
    %                     h(ll) = infCorrectedEntropy(allLatencies,@(x) vasicekEntropy(x,ll,true,2),2);
    %                 end

    %                 figure;
    %                 plot(1:n,h,1:n,infCorrectedEntropy(latencies,@(x) entropy(x',[],'unbiased',2),2)*ones(size(h)));
    % 
    %                 totalLatencyEntropy2(ii,jj,hh,1) = min(h(isfinite(h))); % trololololol
    %                 totalLatencyEntropy2(ii,jj,hh,1) = infCorrectedEntropy(allLatencies,@(x) entropy(x',[],'unbiased',2),2);
    %                 totalLatencyEntropy2(ii,jj,hh,3) = infCorrectedEntropy(1./allLatencies,@(x) entropy(x',[],'unbiased',2),2);

                    if isempty(allValidLatencies)
                        mutualLatencyPhaseInformationPerWidth4(ii,jj,hh,:) = 0;
                        continue;
                    end

                    n = numel(allValidLatencies);

                    % if we don't have enough data points that the probability
                    % of choosing the same (maxSplit)th subset of data is less
                    % than 0.5 (birthday problem), the cell probably isn't
                    % telling you anything about the stimulus anyway
                    if n < maxSplit || nchoosek(n,maxSplit) < repeats^2
                        mutualLatencyPhaseInformationPerWidth4(ii,jj,hh,:) = 0;
                        continue;
                    end

                    nn = 1;
                    for kk = 1:maxSplit
                        m = floor(n/kk);

                        if m < 2 % mutual information is necessarily zero if you only have 0 or 1 samples, so skip these cases
                            mutualLatencyPhaseInformationPerWidth4(ii,jj,hh,nn:end) = nan;
                            break;
                        end

                        for ll = 1:(kk == 1) + (kk > 1)*repeats
                            index = randperm(n);

                            for mm = 1:kk
                                if mm == kk
                                    validLatencies = allValidLatencies(index((mm-1)*m+1:n));
                                    validPhases = allValidPhases(index((mm-1)*m+1:n));
                                else
                                    validLatencies = allValidLatencies(index((mm-1)*m+1:mm*m));
                                    validPhases = allValidPhases(index((mm-1)*m+1:mm*m));
                                end

                                mutualLatencyPhaseInformationPerWidth4(ii,jj,hh,nn) = information(validLatencies',validPhases',[],'unbiased',2);
                                %                 mutualInverseLatencyPhaseInformationPerWidth4(ii,jj,hh) = information(1./allLatencies',allPhases',[],'unbiased',2);

                                nn = nn + 1;
                            end
                        end
                    end
                end
            end
        end

        unbiasedLatencyInfo = zeros(nResponses,nWidths,2);

        figure;
        for ii = 1:nResponses
            clf;

            for jj = 1:nWidths
                for kk = 1:2
                    info1 = squeeze(mutualLatencyPhaseInformationPerWidth4(ii,jj,kk,:));
                    subplot(2,nWidths,nWidths*(kk-1)+jj);
                    hold on;
                    plot(splits,info1,'Color','b','LineStyle','none','Marker','o');
                    info2 = zeros(maxSplit,1);
                    info2(1) = info1(1);

                    ll = 1;
                    for mm = 2:maxSplit
                        info2(mm) = mean(info1(ll+(1:repeats*mm)));
                        ll = ll + repeats*mm;
                    end

                    plot(1:maxSplit,info2,'Color','g');

                    p = polyfit(splits(isfinite(info1)),info1(isfinite(info1)),2);

                    fplot(@(x) polyval(p,x),[0 maxSplit],'Color','r');

                    unbiasedLatencyInfo(ii,jj,kk) = polyval(p,0);
                end
            end
        end
        
        if ee == 1
            unbiasedLatencyInfos = zeros([size(unbiasedLatencyInfo) nRecordings]);
            unbiasedLatencyInfos(:,:,:,1) = unbiasedLatencyInfo;
            indices = 1:size(unbiasedLatencyInfos,1);
        else
            indices = [];
            for ii = 1:nResponses
                if isstruct(cellIndices)
                    indices(end+1) = cellIndices.(cells{ii,1})(cells{ii,2}); %#ok<AGROW>
                else
                    indices(end+1) = cellIndices(cells(ii,1),cells(ii,2),cells(ii,3)+1); %#ok<AGROW>
                end
            end
            
            if max(indices) > size(unbiasedLatencyInfos,1)
                unbiasedLatencyInfos(end+1:max(indices),:,:,:) = 0;
            end
                
            unbiasedLatencyInfos(indices,:,:,ee) = unbiasedLatencyInfo;
        end

    %     mutualLatencyPhaseInformationPerWidth3 = totalLatencyEntropy2 - conditionalLatencyPhaseEntropy2;

        %% MUTUAL INFORMATION BETWEEN LATENCY AND STIMULUS PER PHASE - Moddemeijer method plus updated correction for delta at infinity

    %     mutualLatencyWidthInformationPerPhase = zeros(nResponses,2);
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             allLatencies = [];
    %             allValidLatencies = [];
    %             allValidWidths = [];
    %             
    %             for jj = 1:nWidths
    %                 for kk = 1:nPhases
    %                     latencies = datas{ii,jj,kk,2,hh};
    %                     
    %                     validLatencies = latencies(isfinite(latencies));
    %                     
    %                     allLatencies = [allLatencies; latencies]; %#ok<AGROW>
    %                     allValidLatencies = [allValidLatencies; validLatencies]; %#ok<AGROW>
    %                     allValidWidths = [allValidWidths; widths(jj)*ones(size(validLatencies))]; %#ok<AGROW>
    %                 end
    %             end
    %                 
    %             if isempty(allValidLatencies) || numel(allValidLatencies) == 1
    %                 mutualLatencyWidthInformationPerPhase(ii,hh) = 0;
    %             else
    %                 mutualLatencyWidthInformationPerPhase(ii,hh) = information(allValidLatencies',allValidWidths',[],'unbiased',2);
    %             end
    %         end
    %     end

        %% MUTUAL INFORMATION BETWEEN LATENCY AND STIMULUS - old method
        %* gives nonsense values that cluster around 0.3
    %     
    %     mutualLatencyInformation = zeros(nResponses,nWidths,nPhases,2);
    %     stimulusProbability = 1/(nWidths*nPhases);
    %     bins = 1:maxLatency;
    %     
    %     for hh = 1:2
    %         for ii = 1:nResponses
    %             totalLatencyProbability = zeros(numel(bins),1);
    %             conditionalLatencyProbability = zeros(numel(bins),nWidths,nPhases);
    %             
    %             for jj = 1:nWidths
    %                 for kk = 1:nPhases
    %                     latencies = datas{ii,jj,kk,2,hh};
    %                     [mu,sigma] = normfit(latencies);
    %                     conditionalLatencyProbability(:,jj,kk) = normpdf(bins,mu,sigma);
    %                     totalLatencyProbability = totalLatencyProbability + conditionalLatencyProbability(:,jj,kk)/stimulusProbability;
    %                 end
    %             end
    %             
    %             for jj = 1:nWidths
    %                 for kk = 1:nPhases
    %                     information = stimulusProbability*conditionalLatencyProbability(:,jj,kk).*log(conditionalLatencyProbability(:,jj,kk)./totalLatencyProbability)/log(2);
    %                     information(isnan(information)) = 0; % linear term dominates the multiplication, so 0*-Inf = 0
    %                     mutualLatencyInformation(ii,jj,kk,hh) = mutualLatencyInformation(ii,jj,kk,hh) - trapz(information);
    %                 end
    %             end
    %         end
    %     end

        save([fileDir '\' filename '_info.mat'],'mutualCountPhaseInformationPerWidth2','mutualLatencyPhaseInformationPerWidth4','unbiasedCountInfo','unbiasedLatencyInfo');

        %% MUTUAL INFORMATION PLOTS

    %     fig = figure;
    %     set(fig,'Position',[-1920 -180 1920 1080]);
    %     
    %     mutualCountInformations = {mutualCountWidthInformation mutualCountPhaseInformation mutualCountStimulusInformation};
    %     mutualLatencyInformations = {mutualLatencyWidthInformation mutualLatencyPhaseInformation mutualLatencyStimulusInformation};
    %     titles = {'Bar Width' 'Bar Shift' 'Stimulus Identity'};
    %         
    %     for ii = 1:2
    %         for jj = 1:3
    %             subplot(2,3,3*(ii-1)+jj)
    %             scatter(mutualCountInformations{jj}(:,ii),mutualLatencyInformations{jj}(:,ii))
    %             hold on
    %             fplot(@(x) x,xlim)
    %             xlabel('I{_{count}}/bits');
    %             ylabel('I{_{latency}}/bits');
    %             
    %             if ii == 1
    % %                 title(['O' 'n'*ones(1,ii==1) 'f'*ones(1,2*(ii==2)) ' responses']);
    %                 title(['Mutual information with ' titles{jj}])
    %             end
    %         end
    %     end
    %     
    %     saveas(fig,[fileDir '\' filename '_' opts.latencytype '_latency_info_vs_rate_info.png']);

        fig = figure;
        set(fig,'Position',[-1920 -180 1920 1080]);

        for ii = 1:2
            for jj = 1:nWidths
                subplot(2,nWidths,nWidths*(ii-1)+jj)
    %             scatter(mutualCountPhaseInformationPerWidth(jj,:,ii),mutualLatencyPhaseInformationPerWidth(:,jj,ii));
    %             scatter(mutualCountPhaseInformationPerWidth2(:,jj,ii,2),mutualLatencyPhaseInformationPerWidth4(:,jj,ii));
                scatter(unbiasedCountInfo(:,jj,ii),unbiasedLatencyInfo(:,jj,ii))
                hold on;
                fplot(@(x) x,xlim)

                if ii == 2 && jj == 1
                    xlabel('I{_{count}}/bits');
                    ylabel('I{_{latency}}/bits');
                end

                if ii == 1
                    title(['Width ' num2str(widths(jj)) '{\mu}m']);
                end
            end
        end

    %     saveas(fig,[fileDir '\' filename '_' opts.latencytype '_latency_info_vs_rate_info_per_width.png']);
        saveas(fig,[fileDir '\' filename '_' opts.latencytype '_latency_info_vs_rate_info_per_width_bias_corrected.png']);

    %     for ii = 1:2
    % %         for jj = 1:nPhases
    % %             subplot(2,nPhases,nPhases*(ii-1)+jj)
    %             subplot(1,2,ii);
    % %             scatter(mutualCountPhaseInformationPerWidth(jj,:,ii),mutualLatencyPhaseInformationPerWidth(:,jj,ii));
    %             scatter(mutualCountWidthInformationPerPhase(:,ii,2),mutualLatencyWidthInformationPerPhase(:,ii));
    %             hold on;
    %             fplot(@(x) x,xlim)
    %             
    %             if ii == 11
    %                 xlabel('I{_{count}}/bits');
    %                 ylabel('I{_{latency}}/bits');
    %             end
    %             
    %             if ii == 1
    %                 title(['Phase ' num2str(phase(jj)/nPhases)]);
    %             end
    % %         end
    %     end
    %     
    %     saveas(fig,[fileDir '\' filename '_' opts.latencytype '_latency_info_vs_rate_info_per_phase.png']);
    %     
    %     fig = figure;
    %     set(fig,'Position',[-1920 -180 1920 1080]);
    %     
    %     combinedInformations = {mutualCountStimulusInformation mutualLatencyStimulusInformation};
    %     summedInformations = {mutualCountWidthInformation+mutualCountPhaseInformation mutualLatencyWidthInformation+mutualLatencyPhaseInformation};
    %     variables = {'Count' 'Latency'};
    %     
    %     for ii = 1:2
    %         for jj = 1:2
    %             subplot(2,2,2*(ii-1)+jj);
    %             scatter(summedInformations{ii}(:,jj),combinedInformations{ii}(:,jj));
    %             hold on;
    %             fplot(@(x) x,xlim)
    %             xlabel(['I(' variables{ii} ';Width) + I(' variables{ii} ';Phase)']);
    %             ylabel(['I(' variables{ii} ';(Width,Phase))']);
    %             title(['O' 'n'*ones(1,jj==1) 'f'*ones(1,2*(jj==2)) ' responses']);
    %         end
    %     end
    %     
    %     saveas(fig,[fileDir '\' filename '_combined_vs_summed_info_latencytype_' opts.latencytype '.png']);
    end
    
    [~,currentDir] = fileparts(pwd);
    
    figure;
    hold on;
    
    % TODO : writing indices here only works if there are two recordings
    subplot(1,2,1);
    hold on;
    plot([ones(1,nResponses); 2*ones(1,nResponses)],squeeze(unbiasedCountInfos(indices,1,1,:))','Color','b','Marker','o');
    plot([ones(1,nResponses); 2*ones(1,nResponses)],squeeze(unbiasedCountInfos(indices,2,1,:))','Color','r','Marker','o');
    set(gca,'XTick',[1 2],'XTickLabel',{'No Drug' 'Drug'});
    xlim([0.5 2.5]);
    ylabel('Mutual Information');
    title('Spike Count');
    
    subplot(1,2,2);
    hold on;
    plot([ones(1,nResponses); 2*ones(1,nResponses)],squeeze(unbiasedLatencyInfos(indices,1,1,:))','Color','b','Marker','o');
    plot([ones(1,nResponses); 2*ones(1,nResponses)],squeeze(unbiasedLatencyInfos(indices,2,1,:))','Color','r','Marker','o');
    set(gca,'XTick',[1 2],'XTickLabel',{'No Drug' 'Drug'});
    xlim([0.5 2.5]);
    title('Response Latency');
    
%     plot(1*ones(nCells,1),diff(unbiasedCountInfos(indices,1,1,:),[],4),'Color','b','LineStyle','none','Marker','o');
%     h = plot(3*ones(nCells,1),diff(unbiasedCountInfos(indices,2,1,:),[],4),'Color','b','LineStyle','none','Marker','o');
%     set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     plot(2*ones(nCells,1),diff(unbiasedLatencyInfos(indices,1,1,:),[],4),'Color','r','LineStyle','none','Marker','o');
%     h = plot(4*ones(nCells,1),diff(unbiasedLatencyInfos(indices,2,1,:),[],4),'Color','r','LineStyle','none','Marker','o');
%     set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
%     legend({'Count' 'Latency'});
%     
%     set(gca,'XTick',[1.5 3.5],'XTickLabel',75./widths);
%     
%     xlabel('Bar Speed (mm/s)');
%     ylabel('Mutual Information Change');
%     
%     xlim([0.5 4.5]);
%     
%     line([0.5 4.5],[0 0],'Color','k','LineStyle','--');
%     
    figfile = sprintf('%s_%s',currentDir,opts.infofilesuffix);
%     saveas(gcf,figfile,'fig');
%     saveas(gcf,figfile,'png');
    
    xlsfile = sprintf('%s.xlsx',figfile);
    
    lastRow = 1+nWidths;
    xlswrite(xlsfile,widths(:),'Count',sprintf('A2:A%d',lastRow));
    xlswrite(xlsfile,widths(:),'Latency',sprintf('A2:A%d',lastRow));
    
    nCells = size(unbiasedCountInfos,1);
    header = cell(1,(nRecordings-1)*nCells+2);
    header{1} = 'Width';
    
    for ii = 1:nRecordings
        header((ii-1)*nCells+2) = recordings(ii);
    end
    
    lastColumn = getExcelColumn(size(header,2)-1);
    xlswrite(xlsfile,header,'Count',sprintf('A1:%s1',lastColumn));
    xlswrite(xlsfile,header,'Latency',sprintf('A1:%s1',lastColumn));
    
    unbiasedCountInfos = squeeze(unbiasedCountInfos(:,:,1,:));
    unbiasedCountInfos = permute(unbiasedCountInfos,[2 1 3]);
    unbiasedCountInfos = reshape(unbiasedCountInfos,nWidths,nCells*nRecordings);
    
    lastColumn = getExcelColumn(nCells*nRecordings);
    xlswrite(xlsfile,unbiasedCountInfos,'Count',sprintf('B2:%s%d',lastColumn,lastRow));
    
    unbiasedLatencyInfos = squeeze(unbiasedLatencyInfos(:,:,1,:));
    nCells = size(unbiasedLatencyInfos,1);
    unbiasedLatencyInfos = permute(unbiasedLatencyInfos,[2 1 3]);
    unbiasedLatencyInfos = reshape(unbiasedLatencyInfos,nWidths,nCells*nRecordings);
    
    xlswrite(xlsfile,unbiasedLatencyInfos,'Latency',sprintf('B2:%s%d',lastColumn,lastRow));
end

function Y = pseudolog(X,b)
    if nargin < 2
        b = exp(1);
    end
    
    Y = zeros(size(X));
    Y(X <= 1) = X(X <= 1);
    Y(X > 1) = log(X(X > 1))/log(b) + 1;
end

function h = infCorrectedEntropy(x,entropyFn,base)
    n = numel(x);
    
    xfin = x(isfinite(x));
    
    if isempty(xfin)
        % pinf = 1 so pinf*log(pinf) = 1log 1 = 0 & (1-pinf)*log((1-ping)) = 0log 0 = 0
        h = 0;
        return;
    end
    
    if numel(xfin) == 1
        h = 0;
    else
        h = entropyFn(xfin);
    end
    
    pfin = numel(xfin)/n; 
    
    if pfin == 1
        % pinf*log(pinf) = 0log 0 = 1 & (1-pinf)*log((1-pinf)) = 1log 1 = 0
        return;
    end
    
    % pfin = 1-pinf
    h = pfin*h - (pfin*log(pfin) + (1-pfin)*log(pfin))/log(base);
end

% wrong, don't use
function h = moddemeijerEntropy(x)
    if std(x) == 0
        h = 0;
        return;
    end
    
    xfin = x(isfinite(x));
    xinf = x(isinf(x));
    n = numel(x);
    
    if isempty(xfin) || std(xfin) == 0
        h = 0;
    elseif numel(xfin) == 1
        h = -(1/n)*log2(1/n);
    else
        [f,descriptor] = histogram(xfin');
        f = f/n;
        lower = descriptor(1);
        upper = descriptor(2);
        ncell = descriptor(3);
        
        h = -sum(plogp(f));
        
        h = h+log2(n)+log2((upper-lower)/ncell)+(ncell-1)/(2*n);
    end
    
    if isempty(xinf)
        return;
    end
    
    pinf = numel(xinf)/n;
    h = h - pinf*log2(pinf);
end