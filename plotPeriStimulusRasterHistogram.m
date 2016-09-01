function plotPeriStimulusRasterHistogram(recording,varargin)
    if isstruct(recording) && (~isfield(recording,'spont') || ~all(recording.spont))
        error('calculateSpikeTriggeredAverage only works for continuous MC_Rack data files');
    end

    options = getopt('file=NaN rasterfilesuffix=''abs''',varargin{:});
    
    if ischar(options.file)
        load(options.file,'cells','rasters','histograms','edgess','repeats','stimulusLengths','stimulusTimings','subStimulusLengths','factors','levels','valuess');
        
        [filedir,filename] = fileparts(options.file);
        clustered = true;
    else
        [filedir,filename] = getAnalysisOutputDir(recording);

%         [~,~,~,cells,nCells,~,clustered] = concatenateSpikes(recording,varargin);
        
        load([filedir '\' filename '_psrh_' options.rasterfilesuffix '.mat'],'rasters','histograms','edgess','repeats','stimulusLengths','stimulusTimings','subStimulusLengths','factors','levels','valuess','cells');
        clustered = numel(unique(cells(:,end))) ~= 1; %#ok<NODEF>
    end
    
    nCells = size(cells,1);
    
%     load([filedir '\' recording.dataFile '_vsync_times.mat']);
%     load(recording.stimulusFile,'getExtraParams','getPixels','misses','seed','stimuli','textureRect','timeOffset','vbls','version');
    
    nConditions = prod(levels);
    nFactors = numel(factors); %#ok<USENS>
    % TODO : distinguishable_colors but pale for fills
    colours = 0.1*[1 0 0; 0 1 0; 0 0 1; 1 0 1; 1 1 0; 0 1 1; 0 0 0] + 0.9*ones(7,3);
    titleSuffixes = {' vsyncs before' ' vsyncs after'};
    titleInfixes = {' all spikes ' ' clustered spikes '};
    
    if ~isstruct(recording) || ~isfield(recording,'rasterFn') || ~isa(recording.rasterFn,'function_handle')
        rasterFn = @(varargin) 1;
    else
        rasterFn = recording.rasterFn;
    end
    
    fig(1) = figure('Position',[0 0 1200 900]);
    set(fig(1),'Renderer','zbuffer');
%     fig(2) = figure('Position',[0 0 1200 900]);
    for ii = 1:nCells
        subscript = cell(1,ndims(rasters)); %#ok<USENS>
        subscript{1} = ii;
        
        if iscell(cells)
            channel = cells{ii,1};
            cluster = num2str(cells{ii,2});
        else
            channel = char(cells(ii,1:end-1));
            cluster = num2str(cells(ii,end));
        end
        
        titlePrefix = ['psrh_' options.rasterfilesuffix '_' filename titleInfixes{1+(~isnan(str2double(channel)) && clustered)} 'channel ' channel ' cluster ' cluster];
        
        clf(fig(1));
%         clf(fig(2));

        maxh = -Inf;
        nHistograms = numel(histograms)/nCells; %#ok<USENS>
        
        for jj = 1:nHistograms
            maxh = max(maxh,max(histograms{(jj-1)*nCells+ii}));
        end
        
        tic;
        for jj = 1:nConditions
            [subscript{2:end-1}] = ind2sub(levels,jj);
            
            for kk = 1 %:2
%                 figure(fig(kk));
                
                subscript{end} = kk;
                index = sub2ind(size(rasters),subscript{:});
                
                raster = rasters{index};
                edges = edgess{index}; %#ok<USENS>
                histogram = histograms{index};
                nLines = repeats(index);
                
                [nRows,nCols] = subplots(nConditions);
                rows = ceil(jj/nCols);
                cols = mod(jj-1,nCols)+1;
                position = [(cols-1)/nCols+20/1200 1-rows/nRows+20/900 1/nCols-25/1200 1/nRows-40/900];
                ax = subplot('Position',position);
                set(gca,'FontSize',6);
                hold on;
                
                if exist('stimulusLengths','var')
                    stimulusLengthIndex = sub2ind(size(stimulusLengths),subscript{2:end});
                    maxt = stimulusLengths(stimulusLengthIndex);
                else
                    maxt = 0;
                    subStimulusLengthsSize = size(subStimulusLengths);
                    
                    for ll = 1:subStimulusLengthsSize(end-1)
                        subStimulusLengthSubscript = num2cell([subscript{1+(ndims(subStimulusLengths) == numel(levels)+2):end-1} ll subscript{end}]);
                        subStimulusLengthsIndex = sub2ind(subStimulusLengthsSize,subStimulusLengthSubscript{:});
                        maxt = maxt + subStimulusLengths(subStimulusLengthsIndex);
                    end
                end
                
                if maxh > 0
                    yy = [0 2*maxh];
                else
                    yy = [0 1];
                end
                
                if exist('stimulusTimings','var')
                    timings = stimulusTimings{stimulusLengthIndex}; %#ok<USENS>
                    
                    for ll = 1:size(timings,3)
                        for mm = 1:size(timings,2)
                            fill(timings([1 2 2 1],mm,ll)-timings(1,1,ll),maxh*(1+(ll+[0 0 1 1])/(nLines+1)),colours(mod(mm-1,7)+1,:),'LineStyle','none');
                        end
                    end
                else
                    t = 0;
                    for ll = 1:size(subStimulusLengths,ndims(subStimulusLengths)-1)
                        tt(1) = t;
                        subStimulusLengthSubscript = {subscript{1+(ndims(subStimulusLengths) == numel(levels)+2):end-1} ll subscript{end}};
                        dt = subStimulusLengths(sub2ind(size(subStimulusLengths),subStimulusLengthSubscript{:}));

                        if dt == 0
                            continue;
                        end

                        t = t + dt;
                        tt(2) = t;

                        fill(tt([1 2 2 1]),yy([1 1 2 2]),colours(mod(ll-1,7)+1,:),'LineStyle','none');
                    end
                end
                
                line([0 0],yy,'LineStyle','--','Color',zeros(1,3));
                line([maxt maxt],yy,'LineStyle','--','Color',zeros(1,3));
                
                bar(edges,histogram,'histc');
                
                for ll = 1:nLines
                    if ~isempty(raster{ll})
                        line(repmat(raster{ll}',2,1),maxh*(1+[ll;ll+1]/(nLines+1))*ones(size(raster{ll}))','Color','k');
%                         plot(repmat(raster{ll}',2,1),maxh*(1+[ll;ll+1]/(nLines+1))*ones(size(raster{ll}))','LineStyle','none','Marker','.','Color','k');
                    end
                end
                
                rasterFn(ax,yy,recording,channel,cluster,subscript,valuess); %#ok<USENS>
                
                xlim([-1 maxt+1]);
                ylim(yy);
                
                xlabel('Time/s');
                ylabel('Firing Rate/Hz');
                plotTitle = '';
                
                for ll = 1:nFactors
                    values = valuess{ll};
                    value = values(subscript{1+ll},:);
                    factor = factors{ll};
                    
                    if ~iscell(factor)
                        subfactors = {factor};
                    else
                        subfactors = factor;
                    end
                        
                    for mm = 1:numel(subfactors)
                        subfactor = subfactors{mm};
                        
                        if numel(value) == numel(subfactors)
                            subvalue = value(mm);
                        else
                            subvalue = value;
                        end

                        if isnumeric(subvalue)
                            subvalue = num2str(horzcat(subvalue));
                        end
                        
                        plotTitle = [plotTitle ' ' subfactor ' ' subvalue]; %#ok<AGROW>
                    end
                end
                
                title(plotTitle,'FontSize',6);
            end
        end
        toc;
        
        for kk = 1 %:2
            figure(fig(kk));
            set(fig(kk), 'PaperPositionMode', 'auto')
            figureTitle = [titlePrefix titleSuffixes{kk}];
%             mtit(figureTitle);
            figfile = [filedir '\' strrep(figureTitle,' ','_')];
%             tic;
%             saveas(fig(kk),[filename '.fig']);
%             toc;
            tic;
            saveas(fig(kk),[figfile '.png']);
            toc;
        end
    end
end