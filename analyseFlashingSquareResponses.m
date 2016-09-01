function analyseFlashingSquareResponses(stimRecording,spontRecording,varargin)
    [fileDir,filename] = getAnalysisOutputDir(stimRecording);
    filePrefix = [fileDir '\' filename];
    load([filePrefix '_psrh_abs.mat']);
    
    if numel(factors) ~= 1 || numel(valuess) ~= 1 || ~strcmp(factors{1},'rect') %#ok<USENS>
        error('File %i (%s) was not a flashing spots experiment',filename);
    end
    
    nCells = size(repeats,1);
    
    % TODO : hack.  is this even needed any more?
    try
        spontDir = getAnalysisOutputDir(spontRecording);
        load([spontDir '\SNR.mat']);
    catch
        load([fileDir '\SNR.mat']);
    end
    
    T = 0;
    nStimuli = 0;
    nTrials = -Inf;
    for ii = 1:numel(stimulusTimings) %#ok<USENS>
        timings = stimulusTimings{ii};
        n = size(timings,3);
        nTrials = max(nTrials,n);
        nStimuli = nStimuli + n;
        
        for jj = 1:size(timings,3)
            T = T + timings(2,2,jj) - timings(1,1,jj);
        end
    end
    
    approximateStimulusLength = T/nStimuli;
    
    load([filePrefix '_photodiode_timings.mat'],'recordingStartTime');
    [sortedSpikeTimes,sortedSpikeClusters,sortedSpikeChannels,cells] = concatenateSpikes(stimRecording,varargin);
    
    cellTypesFile = sprintf('%s\\%s_cell_types.mat',fileDir,filename);
    if exist(cellTypesFile,'file')
        load(cellTypesFile,'cellTypes');
        responsive = isfinite(cellTypes(:,3)); %#ok<NODEF>
        responsiveTypes = find(responsive);
        responsiveCells = find(ismember([str2num(char(cells(:,1:2))) cells(:,3)],cellTypes(responsive,1:2),'rows'));
    else
        responsiveCells = 1:nCells;
    end
    
    nResponses = numel(responsiveCells);
    
    nonstationary = zeros(nResponses,2);
    
    figure;
    for ii = 1:nResponses
        index = responsiveCells(ii);
        spikeTimes = sortedSpikeTimes(sortedSpikeClusters == cells(index,3) & sortedSpikeChannels(:,1) == cells(index,1) & sortedSpikeChannels(:,2) == cells(index,2));
        
        nSpikes = zeros(nStimuli,1);
        
        for jj = 1:nStimuli
            nSpikes(jj) = sum(spikeTimes >= (jj-1)*approximateStimulusLength & spikeTimes < jj*approximateStimulusLength);
        end
        
%         [~,nonstationary(ii,1)] = kstest2(nSpikes(1:nStimuli/2),nSpikes(nStimuli/2+1:end));
%         
%         meanNSpikes = mean(reshape(nSpikes,nTrials,nStimuli/nTrials))';
%         
%         [~,nonstationary(ii,2)] = kstest2(meanNSpikes(1:nStimuli/(2*nTrials)),meanNSpikes(nStimuli/(2*nTrials)+1:end));
        plot(0:approximateStimulusLength:T-approximateStimulusLength,nSpikes)
    end
    
    rects = valuess{1};
    
    minX = min(rects(:,1));
    minY = min(rects(:,2));
    
    width = diff(rects(1,[1 3]));
    
    maxX = max(rects(:,1))+width;
    maxY = max(rects(:,2))+width;
    
    assert( ...
        mean(rects(:,3)-rects(:,1)) == width && ...
        mean(rects(:,4)-rects(:,2)) == width && ...
        std(rects(:,3)-rects(:,1)) == 0 && ...
        std(rects(:,4)-rects(:,2)) == 0, ...
        'All square sizes should be the same' ...
        );
    
    pixels = zeros(maxY-minY,maxX-minX,nResponses);
    
    load(['./' filename '.mat']);
    nTiles = sqrt(size(getExtraParams,2));
    
    tileWidth = maxX-minX;
    
    tileGrid = tileWidth*[repmat(0:nTiles-1,1,nTiles)' kron(0:nTiles-1,ones(1,nTiles))'];
    
    opts = getopt('locsonly=''no'' ratetype=''abs'' responselen=0.5 square=NaN',varargin{:});
    square = opts.square;
    if isnumeric(square) && numel(square) == 4 && all(isfinite(square))
        assert(square(3) == square(4),'That''s not a square');
    else
        square = [minX minY 60 60];
    end
    
    micronsPerPixel = 1430/square(3);
    electrodeRadius = ((30/micronsPerPixel)/2);
    electrodeSpacing = (200/micronsPerPixel);
    electrodeCoords = electrodeRadius+(0:7)*electrodeSpacing;
    electrodeXs = square(1)+electrodeCoords;
    electrodeYs = square(2)+electrodeCoords;
    
%     load(sprintf('%s\\%s_isi_stats',fileDir,filename),'lambdas');
    
    responseLocations = cell(60,max(cells(:,3)),2); %#ok<NODEF>
    responseStrengths = cell(60,max(cells(:,3)),2);
%     responsiveCells = [];
%     responsePolarities = [];
    
    channels = zeros(nResponses,1);
    clusters = zeros(nResponses,1);
    coords = nan(nResponses,2);
%     polarities = nan(nResponses,1);

    k = [-1 1; 1 1; 1 -1];
    
    for ii = 1:nResponses
        index = responsiveCells(ii);
        channel = cells(index,1:2);
        channelLabel = str2double(char(channel));
        channelIndex = mcsChannelNumberToChannelIndex(channelLabel);
        
        channels(ii) = channelLabel;
        
        cluster = cells(index,3);
        
        clusters(ii) = cluster;
        
        if cluster == 0
            continue;
        end
        
        allNSpikes = [];
        meanNSpikes = zeros(levels,2);
        
        polarity = cellTypes(responsiveTypes(ii),3);
        
        nn = 1;
        for jj = 1:levels
            stimulusTimes = stimulusTimings{jj,1};
            raster = rasters{index,jj,1}; %#ok<USENS>
            nTrials = numel(raster);
            
            nSpikess = zeros(nTrials,2);
            stimLengths = zeros(nTrials,2);
            
            for kk = 1:nTrials
                for ll = 1:2
                    interval = stimulusTimes(1,ll,kk)-stimulusTimes(1,1,kk)+[0 opts.responselen];
                    
                    nSpikes = sum(raster{kk} >= interval(1) & raster{kk} < interval(2));
                    stimLength = diff(interval);
                    
                    if strcmp(opts.ratetype,'abs') || strcmp(opts.ratetype,'rel')
                        evokedRate = nSpikes/stimLength;
                    
                        if strcmp(opts.ratetype,'rel')
                            if isempty(spikeRates{channelIndex}) || cluster > numel(spikeRates{channelIndex}) %#ok<USENS>
                                spontRate = 0;
                            else
                                spontRate = spikeRates{channelIndex}(cluster);
                            end

                            if spontRate == 0
                                rateChange = evokedRate;
                            else
                                rateChange = (evokedRate-spontRate)/spontRate;
                            end
                        else
                            rateChange = evokedRate;
                        end
                    
                        meanNSpikes(jj,ll) = meanNSpikes(jj,ll) + rateChange/nTrials;
                    elseif strcmp(opts.ratetype,'delta')
                        preSpikes = sum( ...
                            sortedSpikeTimes+recordingStartTime >= stimulusTimes(1,ll,kk)-opts.responselen & ...
                            sortedSpikeTimes+recordingStartTime < stimulusTimes(1,ll,kk) & ...
                            sortedSpikeChannels(:,1) == channel(1) & ...
                            sortedSpikeChannels(:,2) == channel(2) & ...
                            sortedSpikeClusters == cluster ...
                            );
                        
                        if nSpikes == 0 && preSpikes == 0
                            evokedRate = 1;
                        else
                            evokedRate = nSpikes/preSpikes-1;
                        end 
                        
                        rateChange = evokedRate;
                    else
                        error('Unknown rate type %s\n',opts.ratetype);
                    end
                    
                    xs = rects(jj,1)-minX+(1:width);
                    ys = rects(jj,2)-minY+(1:width);
                    
                    pixels(xs,ys,ii) = pixels(xs,ys,ii) + k(polarity+2,ll)*evokedRate/nTrials; 
%                     pixels(xs,ys,ii,ll) = pixels(xs,ys,ii,ll) + rateChange/nTrials; 
                    nSpikess(kk,ll) = rateChange; %nSpikes;   
                    stimLengths(kk,ll) = stimLength;
                    
                    allNSpikes(end+1) = rateChange; %#ok<AGROW>
                end
            end
            
%             try
%                 lambda = lambdas{channelIndex}(cluster); %#ok<USENS>
%             catch err %#ok<NASGU>
%                 lambda = 0;
%             end
            if strcmp(opts.ratetype,'delta')
                meanNSpikes(jj,:) = median(nSpikess);
            end
            
            nSpikes = median(nSpikess);
            stimLength = mean(stimLengths);

%             for ll = 1:2
%                 if (lambda == 0 && nSpikes(ll) > 0) || (lambda > 0 && (1-poisscdf(ceil(nSpikes(ll)-1),lambda*stimLength(ll))) < 0.05)
%                     xy = rects(jj,1:2)-[minX minY]+width/2;
%                     
%                     if isempty(responseLocations{channelIndex,cluster,ll})
%                         responseLocations{channelIndex,cluster,ll} = xy;
%                         responseStrengths{channelIndex,cluster,ll} = nSpikes(ll);
%                     else
%                         responseLocations{channelIndex,cluster,ll}(end+1,:) = xy;
%                         responseStrengths{channelIndex,cluster,ll}(end+1) = nSpikes(ll);
%                     end
%                 end
%             end
        end
        
        if strcmp(opts.locsonly,'yes')
            [~,maxI] = max(reshape(meanNSpikes,levels*2,1));

            [location,polarity] = ind2sub([levels 2],maxI);

            polarities(ii) = 3-2*polarity;

            possibleCoords = repmat(rects(location,1:2),nTiles^2,1)+tileGrid+width/2;

            % need to compensate for rotation and reflection of screen: see
            % lab book #2, experiment JBWT0045 for details.  Also note that
            % the y-axis will be flipped in the plots.
            electrodeX = electrodeXs(9-mod(channelLabel,10));
            electrodeY = electrodeYs(9-floor(channelLabel/10));

            distances = sqrt(sum((possibleCoords-repmat([electrodeX electrodeY],nTiles^2,1)).^2,2));

            [~,minI] = min(distances);
            coords(ii,:) = possibleCoords(minI,:);

            continue;
        end
        
%         sample = bootstrap(allNSpikes,10,10000,true);
%         ps = zeros([size(meanNSpikes) 2]);
%         
%         for jj = 1:levels
%             for kk = 1:2
%                 ps(jj,kk,1) = sum(sample >= meanNSpikes(jj,kk))/numel(sample);
%                 ps(jj,kk,2) = sum(sample <= meanNSpikes(jj,kk))/numel(sample);
%             end
%         end
%         
%         alpha = 0.05/(nResponses*2); %*levels*2);
%         
%         response = squeeze(any(ps <= alpha));
%         
%         if any(any(response))
%             responsiveCells(end+1) = ii; %#ok<AGROW>
%             responsePolarities(end+1) = ...
%                 (response(1,1) || response(2,2)) - ...
%                 (response(2,1) || response(1,2)); %#ok<AGROW>
%         end
    end
    
    if strcmp(opts.locsonly,'yes')
        save(sprintf('%s\\%s_cell_types.mat',fileDir,filename),'channels','clusters','coords','polarities');

        figure;
        scatter(coords(:,1),coords(:,2));
        hold on;

        for ii = 1:numel(electrodeXs)
            for jj = 1:numel(electrodeYs)
                if all(ismember([ii jj],[1 numel(electrodeCoords)]))
                    continue;
                end

                x0 = electrodeXs(ii);
                y0 = electrodeYs(jj);
                r = electrodeRadius;

                [x,y] = ellipseFn(x0,y0,r,r,0);
                h = ezplot(x,y);
                set(h,'Color','k');
            end
        end

        text(electrodeXs(1),  electrodeYs(1),  '88','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(electrodeXs(1),  electrodeYs(end),'18','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(electrodeXs(end),electrodeYs(1),  '81','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(electrodeXs(end),electrodeYs(end),'11','FontSize',6,'HorizontalAlignment','center','VerticalAlignment','middle');

        xlim([0 640]);
        ylim([0 480]);

        figfile = sprintf('%s\\%s_cell_locations',fileDir,filename);
        saveas(gcf,figfile,'fig');
        saveas(gcf,figfile,'png');

        return;
    end
    
%     nResponses = numel(responsiveCells);
%     pixels = pixels(:,:,responsiveCells,:);
    pixels(isinf(pixels)) = max(pixels(isfinite(pixels)));
    tileWidth = size(pixels,2)/width;

    
    
    stas = nan(size(pixels).*[nTiles nTiles 1]./[width width 1]);
    rfs = zeros(size(pixels).*[nTiles nTiles 1]);
    
    rfParams = struct( ...
        'centreX',        cell(nResponses,1), ...
        'centreY',        cell(nResponses,1), ...
        'sdCentreX',      cell(nResponses,1), ...
        'sdCentreY',      cell(nResponses,1), ...
        'sdSurroundX',    cell(nResponses,1), ...
        'sdSurroundY',    cell(nResponses,1), ...
        'theta',          cell(nResponses,1), ...
        'csRatio',        cell(nResponses,1), ...
        'polarity',         cell(nResponses,1)  ...
        );
%         'centreXOn',        cell(nResponses,1), ...
%         'centreXOff',       cell(nResponses,1), ...
%         'centreYOn',        cell(nResponses,1), ...
%         'centreYOff',       cell(nResponses,1), ...
%         'sdCentreXOn',      cell(nResponses,1), ...
%         'sdCentreYOn',      cell(nResponses,1), ...
%         'sdSurroundXOn',    cell(nResponses,1), ...
%         'sdSurroundYOn',    cell(nResponses,1), ...
%         'sdCentreXOff',     cell(nResponses,1), ...
%         'sdCentreYOff',     cell(nResponses,1), ...
%         'sdSurroundXOff',   cell(nResponses,1), ...
%         'sdSurroundYOff',   cell(nResponses,1), ...
%         'thetaOn',          cell(nResponses,1), ...
%         'thetaOff',         cell(nResponses,1), ...
%         'csRatioOn',        cell(nResponses,1), ...
%         'csRatioOff',       cell(nResponses,1), ...
%         'polarity',         cell(nResponses,1)  ...
%         );
    
    electrodeSpacing = electrodeSpacing/width;
    electrodeRadius = electrodeRadius/width;
    arrayWidth = 60/width;
    arrayEdge = (size(stas,1)-arrayWidth)/2;
    electrodeCoords = arrayEdge+electrodeRadius+(0:7)*electrodeSpacing;
    
    avgRFDiameter = 300;
    avgSD = avgRFDiameter/(micronsPerPixel*width*2*2);
    
    maxRFDiameter = 1000;
    maxSD = maxRFDiameter/(micronsPerPixel*width*2*2);
    
    fminconOptions = optimset('Algorithm','sqp','MaxFunEvals',1000);
    
    staFigure = figure;
    rfFigure = figure;
    suffixes = {'On' 'Off'};
    for ii = 1:nResponses
%         polarity = responsePolarities(ii);
        polarity = cellTypes(responsiveTypes(ii),3);
        
%         switch polarity
%             case -1
%                 extremum = @min;
%             otherwise
%                 extremum = @max;
%         end
        extremum = @max;        
        
%         for jj = 1:2
%             if (jj == 1 && polarity == -1) || (jj == 2 && polarity == 1)
%                 continue;
%             end
            
            pix = repmat(pixels(:,:,ii),nTiles,nTiles);
            pix = pix(1:width:end,1:width:end);
            
            figure(staFigure);
            clf;
            hold on;
            surf(pix);
            
            try
                caxis([min(min(pix)) max(max(pix))]);
            catch err %#ok<NASGU>
                caxis([0 1]);
            end
            
            colormap(gray);
            view(2);
            shading flat;
            
            chlabel = cells(responsiveCells(ii),1:2)-48;
            
            % need to compensate for rotation and reflection of screen: see
            % lab book #2, experiment JBWT0045 for details.  Also note that
            % the y-axis will be flipped in the plots.
            electrodeX = electrodeCoords(9-chlabel(2));
            electrodeY = electrodeCoords(9-chlabel(1));
            
            [x,y] = ellipseFn(electrodeX,electrodeY,electrodeRadius,electrodeRadius);
            h = ezplot(x,y);
            set(h,'Color','m','ZData',2*max(max(pix))*ones(size(get(h,'XData'))));
            
            xlim([1 size(pix,2)]);
            ylim([1 size(pix,1)]);
            
            maxI = find(pix == extremum(extremum(pix)));
            
            minDistance = Inf;
            
            oldX0 = NaN;
            oldY0 = NaN;
            
            for kk = 1:numel(maxI)
                [y0,x0] = ind2sub(size(pix),maxI(kk));
                
                distance = sqrt((electrodeX-x0).^2+(electrodeY-y0).^2);
                
                if distance < minDistance
                    minDistance = distance;
                    oldX0 = x0;
                    oldY0 = y0;
                end
            end
            
            x0 = oldX0;
            y0 = oldY0;
            
            if nTiles > 1
                xs = max(1,floor(electrodeX-tileWidth/2)):min(size(pix,1),ceil(electrodeX+tileWidth/2));
                ys = max(1,floor(electrodeY-tileWidth/2)):min(size(pix,1),ceil(electrodeY+tileWidth/2));

                line(xs([1 end end 1; end end 1 1]),ys([1 1 end end; 1 end end 1]),2*max(max(pix))*ones(2,4),'Color','g');
            
                while x0 < tileWidth/2+1
                    x0 = oldX0+tileWidth;
                end

                while x0+tileWidth/2 > size(pix,2)
                    x0 = oldX0-tileWidth;
                end

                while y0 < tileWidth/2+1
                    y0 = oldY0+tileWidth;
                end

                while y0+tileWidth/2 > size(pix,2)
                    y0 = oldY0-tileWidth;
                end

                xs = x0+(-tileWidth/2:tileWidth/2);
                ys = y0+(-tileWidth/2:tileWidth/2);

                line(xs([1 end end 1; end end 1 1]),ys([1 1 end end; 1 end end 1]),2*max(max(pix))*ones(2,4),'Color','b');
                line(x0+[0 1 1 0; 1 1 0 0],y0+[0 0 1 1; 0 1 1 0],2*max(max(pix))*ones(2,4),'Color','y');
            else
                xs = 1:size(pix,2);
                ys = 1:size(pix,1);
            end
    
            figfile = sprintf('%s\\%s_channel_%d%d_cluster_%d_flashing_square_sta',fileDir,filename,chlabel(1),chlabel(2),cells(responsiveCells(ii),3));
            saveas(staFigure,figfile,'fig');
            saveas(staFigure,figfile,'png');

            rfPolarity = (max(max(pix)) > 0) - (min(min(pix)) < 0);
%             if rfPolarity == 0;
%                 pix2 = zeros(size(pix)); %median(reshape(pix,numel(pix),1))*ones(size(pix));
%             else
                pix2 = zeros(size(pix));
%             end

%             pix2 = zeros(size(pix));
            pix2(ys,xs) = pix(ys,xs)-median(reshape(pix,numel(pix),1));
            
            figure(rfFigure);
            clf;
            hold on;
            h = surf(pix2);
            
            try
                caxis([min(min(pix2)) max(max(pix2))]);
            catch err %#ok<NASGU>
                caxis([0 1]);
            end
            
            view(3);
            
            xlim([1 size(pix2,2)]);
            ylim([1 size(pix2,1)]);
            
            if max(max(pix2)) > min(min(pix2))
                zlim([min(min(pix2)) max(max(pix2))]);
            end
            
            [Y,X] = ndgrid(1:size(pix2,1),1:size(pix2,2));
            [Yfine,Xfine] = ndgrid(1/width:1/width:size(pix2,1),1:size(pix2,2));
            
            % TODO : need a way to force same eccentricity in centre and
            % surround
            if false && rfPolarity == 0;
                rfFn = @(b) dog(X-x0-0.5,Y-y0-0.5,b(1),b(2),b(3),NaN,b(4),b(5),true);
                
                b0 = [avgSD*[1 1] 2 0 0.5];
%                 b0 = [avgSD*[1 1 2 2] 0 0.5];
                A = []; %A = [1 0 -1 0 0 0; 0 1 0 -1 0 0; -2 0 1 0 0 0; 0 -2 0 1 0 0];
                B = []; %B = zeros(4,1);
                bmin = [repmat(1/width,1,2) 1 0 0];
                bmax = [repmat(maxSD,1,2) 5 2*pi 1];
            else
                rfFn = @(b) gauss2d(X-x0,Y-y0,b(1),b(2),b(3),true);
                
                b0 = [repmat(avgSD,1,2) 0];
                A = [1 -2 0; -2 1 0];
                B = [0; 0];
                bmin = [repmat(eps,1,2) 0];
                bmax = [repmat(maxSD,1,2) 2*pi];
            end
            
%             if polarity == -1
%                 k = -1;
%             else
%                 k = 1;
%             end
            k = 1;

            mse = @(b) mean(reshape((pix2-k*extremum(extremum(pix2))*rfFn(b)).^2,numel(pix2),1));
            
            b = fmincon(mse,b0,A,B,[],[],bmin,bmax,[],fminconOptions);
            
%             rf = k*extremum(extremum(pix2))*rfFn(b);
%             rf = k*extremum(extremum(pix2))*dog(Xfine-x0-0.5,Yfine-y0-0.5,b(1),b(2),b(3),NaN,b(4),b(5),true);
            rf = k*extremum(extremum(pix2))*gauss2d(Xfine-x0,Yfine-y0,b(1),b(2),b(3),true);
            
            surf(Xfine,Yfine,rf);
            shading interp;
            set(h,'EdgeColor','k','FaceColor','none');
            
            figfile = sprintf('%s\\%s_channel_%d%d_cluster_%d_flashing_square_rf',fileDir,filename,chlabel(1),chlabel(2),cells(responsiveCells(ii),3));
            saveas(rfFigure,figfile,'fig');
            saveas(rfFigure,figfile,'png');
            
            Xoffset = oldX0-x0;
            xs = xs+Xoffset;
            xs = xs(xs >= 1 & xs <= size(pix2,2));
            
            Yoffset = oldY0-y0;
            ys = ys+Yoffset;
            ys = ys(ys >= 1 & ys <= size(pix2,1));
            
            stas(xs,ys,ii) = pix(xs,ys);
            
            x0 = (oldX0+0.5)*width;
            y0 = (oldY0+0.5)*width;
            
%             b = b.*[repmat(width,1,2) 1 1 1];
            b = b.*[repmat(width,1,2) 1];
            
            [X,Y] = ndgrid((1:size(rfs,1))-y0,(1:size(rfs,2))-x0);
%             rfs(:,:,ii) = dog(X,Y,b(1),b(2),b(3),NaN,b(4),b(5),true);
            rfs(:,:,ii) = gauss2d(X,Y,b(1),b(2),b(3),true);
            
            suffix = '';
            rfParams(ii).(sprintf('centreX%s',suffix))      = x0;
            rfParams(ii).(sprintf('centreY%s',suffix))      = y0;
            rfParams(ii).(sprintf('sdCentreX%s',suffix))    = b(1);
            rfParams(ii).(sprintf('sdCentreY%s',suffix))    = b(2);
%             rfParams(ii).(sprintf('sdSurroundX%s',suffix))  = b(3)*b(1);
%             rfParams(ii).(sprintf('sdSurroundY%s',suffix))  = b(3)*b(2);
            rfParams(ii).(sprintf('theta%s',suffix))        = b(3);
%             rfParams(ii).(sprintf('k%s',suffix))            = b(5);
%             rfParams(ii).polarity = polarity;
            rfParams(ii).polarity = cellTypes(responsiveTypes(ii),3);
%         end
    end
    
    responsiveChannels = str2double(cellstr(char(cells(responsiveCells,1:2)))); %#ok<NASGU>
    responsiveClusters = cells(responsiveCells,3); %#ok<NASGU>
    save(sprintf('%s\\%s_flashing_squares',fileDir,filename),'stas','rfs','rfParams','responsiveChannels','responsiveClusters','responsiveCells');
    
    figure;
    hold on;
    
    for ii = 1:numel(electrodeCoords)
        for jj = 1:numel(electrodeCoords)
            if all(ismember([ii jj],[1 numel(electrodeCoords)]))
                continue;
            end
            
            x0 = electrodeCoords(ii)*width;
            y0 = electrodeCoords(jj)*width;
            r = electrodeRadius*width;
            
            [x,y] = ellipseFn(x0,y0,r,r,0);
            h = ezplot(x,y);
            set(h,'Color','k');
        end
    end
    
    text(width*electrodeCoords(1),  width*electrodeCoords(1),  '88','HorizontalAlignment','center','VerticalAlignment','middle');
    text(width*electrodeCoords(1),  width*electrodeCoords(end),'18','HorizontalAlignment','center','VerticalAlignment','middle');
    text(width*electrodeCoords(end),width*electrodeCoords(1),  '81','HorizontalAlignment','center','VerticalAlignment','middle');
    text(width*electrodeCoords(end),width*electrodeCoords(end),'11','HorizontalAlignment','center','VerticalAlignment','middle');
        
    minX = (electrodeCoords(1)-electrodeRadius)*width;
    minY = minX;
    maxX = (electrodeCoords(end)+electrodeRadius)*width;
    maxY = maxX;
    
    for ii = 1:nResponses
        params = rfParams(ii);
        
%         if params.polarity == -1
%             suffixes = {'Off'};
%         elseif params.polarity == 0
%             suffixes = {'On' 'Off'};
%         elseif params.polarity == 1
%             suffixes = {'On'};
%         else
%             continue;
%         end
        suffixes = {''};
        
        for jj = 1:numel(suffixes)
            suffix = suffixes{jj};
            x0 = params.(sprintf('centreX%s',suffix));
            y0 = params.(sprintf('centreY%s',suffix));
            
            Ac = 2*params.(sprintf('sdCentreX%s',suffix));
            Bc = 2*params.(sprintf('sdCentreY%s',suffix));
%             As = 2*params.(sprintf('sdSurroundX%s',suffix));
%             Bs = 2*params.(sprintf('sdSurroundY%s',suffix));
            
%             minX = min(minX,min(x0-As,x0-Bs));
%             minY = min(minY,min(y0-As,x0-Bs));
%             
%             maxX = max(maxX,max(x0+As,x0+Bs));
%             maxY = max(maxY,max(y0+As,x0+Bs));
            
            minX = min(minX,min(x0-Ac,x0-Bc));
            minY = min(minY,min(y0-Ac,x0-Bc));
            
            maxX = max(maxX,max(x0+Ac,x0+Bc));
            maxY = max(maxY,max(y0+Ac,x0+Bc));
            
            theta = rfParams(ii).(sprintf('theta%s',suffix));
            
            [x,y] = ellipseFn(x0,y0,Ac,Bc,theta);
            h = ezplot(x,y);
            set(h,'Color','g');
            
%             [x,y] = ellipseFn(x0,y0,As,Bs,theta);
%             h = ezplot(x,y);
%             set(h,'Color','r');
        end
    end
    
    xlim([minX maxX]);
    ylim([minY maxY]);
    
    figfile = sprintf('%s\\%s_rf_locations',fileDir,filename);
    saveas(gcf,figfile,'fig');
    saveas(gcf,figfile,'png');
    
    close all;
end