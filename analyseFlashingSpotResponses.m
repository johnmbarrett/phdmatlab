function analyseFlashingSpotResponses(recording,varargin)
    fileDir = getAnalysisOutputDir(recording);
    filePrefix = [fileDir '\' recording.dataFile];
    load([filePrefix '_psrh.mat']);
    
    if numel(factors) ~= 1 || numel(valuess) ~= 1 || ~strcmp(factors{1},'rect') %#ok<USENS>
        error('File %i (%s) was not a flashing spots experiment',recording.index,recording.dataFile);
    end
    
    load([filePrefix '_photodiode_timings.mat']);
    [sortedSpikeTimes,sortedSpikeClusters,sortedSpikeChannels,cells,nCells] = concatenateSpikes(recording,varargin);
    
    maxNClusters = max(sortedSpikeClusters)+1;
    figures = zeros(maxNClusters,4);
    dirtyFigures = false(size(figures));
    
    for ii = 1:maxNClusters*4
        figures(ii) = figure;
    end
    
    cleanupFigures = onCleanup(@() close(figures));
    
    if ~isfield(recording,'textureRect')
        load([recording.stimulusFile '.mat'],'textureRect');
    else
        textureRect = recording.textureRect;
    end
    
    X = textureRect(1);
    Y = textureRect(2);
    W = diff(textureRect([1 3]));
    H = diff(textureRect([1 3]));
    
    colours = 'rb';
    
    for ii = 1:nCells
        channel = cells(ii,1:2);
        cluster = cells(ii,3);
        
        baselineEnd = floor((stimulusTimes(1)-recordingStartTime)*10)/10;
        
        edges = linspace(0,baselineEnd,baselineEnd*10+1);
        
        spikeTimes = sortedSpikeTimes( ...
            sortedSpikeTimes < baselineEnd & ...
            strcmp(cellstr(char(sortedSpikeChannels)),char(channel)) & ...
            sortedSpikeClusters == cluster);
        
        firingRate = histc(spikeTimes,edges);
        
        mu = mean(firingRate);
        sigma = std(firingRate);
        
        responses = {[] []};
        peakFRs = {[] []};
        peakLatencies = {[] []};
        polarities = {[] []};
        
        for jj = 1:levels
            edges = edgess{ii,jj,1}; %#ok<USENS>
            minT = find(edges > 0,1);
            maxT = find(edges <= stimulusLengths(jj,1),1,'last');
            offT = median(diff(squeeze(stimulusTimings{jj,1}(1,[1 2],:)))); %#ok<USENS>
            
            histogram = histograms{ii,jj,1}(minT:maxT); %#ok<USENS>
            
            if max(histogram) > mu+2.5*sigma
                responses{1} = [responses{1}; jj];
                peakFRs{1} = [peakFRs{1}; max(histogram)];
                maxLatency = edges(find(histogram == max(histogram),1));
                peakLatencies{1} = [peakLatencies{1}; maxLatency];
                polarities{1} = [polarities{1}; maxLatency < offT];
            end
            
            if min(histogram) < mu-2.5*sigma
                responses{2} = [responses{2}; jj];
                peakFRs{2} = [peakFRs{2}; min(histogram)];
                minLatency = edges(find(histogram == min(histogram),1));
                peakLatencies{2} = [peakLatencies{2}; minLatency];
                polarities{2} = [polarities{2}; minLatency < offT];
            end
        end
        
        if isempty(responses{1}) && isempty(responses{2})
            continue;
        end

        pixels = zeros(H,W);
        subplotnum = str2double(char(channel(1)))+8*str2double(char(channel(2)))-8;
        
        for jj = 1:2
            if isempty(responses{jj})
                continue;
            end
            
            responsiveRects = valuess{1}(responses{jj},:);
            
            responsiveCircles = [mean(responsiveRects(:,[1 3]),2) mean(responsiveRects(:,[2 4]),2) mean([diff(responsiveRects(:,[1 3]),[],2) diff(responsiveRects(:,[1 3]),[],2)],2)];
            responsiveCircles = responsiveCircles - repmat([X Y 0],size(responsiveCircles,1),1);
            
            [sortedRadii,sortedIndices] = sort(responsiveCircles(:,3));
            sortedFRs = peakFRs{jj}(sortedIndices);
            sortedLatencies = peakLatencies{jj}(sortedIndices);
            
            figure(figures(cluster+1,jj));
            subplot(8,8,subplotnum);
            [~,h1,h2] = plotyy(sortedRadii,sortedFRs,sortedRadii,sortedLatencies);
            set(h1,'LineStyle','none','Marker','o');
            set(h2,'LineStyle','none','Marker','o');
            dirtyFigures(cluster+1,jj) = true;

            figure(figures(cluster+1,3));
            subplot(8,8,subplotnum);
            hold on;

            for kk = 1:size(responsiveCircles,1)
                x0 = responsiveCircles(kk,1);
                y0 = responsiveCircles(kk,2);
                r = responsiveCircles(kk,3);
                fplot(@(x) y0 + sqrt(y0^2-(x.^2-2*x0*x+x0^2+y0^2-r^2))*[-1 1],[x0-r x0+r],'Color',colours(jj))

                [xs,ys] = meshgrid(1:W,1:H);
                spot = (xs-x0).^2 + (ys-y0).^2 <= r.^2;

                pixels = pixels+(2*polarities{jj}(kk)-1)*(peakFRs{jj}(kk)-mu)*spot';
            end

            xlim([0 W]);
            ylim([0 H]);
            
            dirtyFigures(cluster+1,3) = true;
        end
        
        maxP = max(max(pixels));
        minP = min(min(pixels));
        
        pixels = 255*(pixels-minP)/(maxP-minP);
        
        figure(figures(cluster+1,4));
        subplot(8,8,subplotnum);
        image(pixels);
        colormap(gray(256));
        
        dirtyFigures(cluster+1,4) = true;
    end
    
    plotnames = {'onpeakfrlat' 'offpeakfrlat' 'stimuli' 'receptive_fields'};
    
    for ii = 1:maxNClusters
        for jj = 1:4
            if ~dirtyFigures(ii,jj)
                continue;
            end
            
            set(figures(ii,jj),'Position',[0 0 4096 4096]);
            filename = [filePrefix '_spotanal_' plotnames{jj} '_cluster_' num2str(ii-1)];
            saveas(figures(ii,jj),filename,'fig');
            saveas(figures(ii,jj),filename,'png');
        end
    end
end