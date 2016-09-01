function rasterHistogramPlot(varargin)
    if nargin == 4
        ax = varargin{1};
        raster = varargin{2};
        histogram = varargin{3};
        edges = varargin{4};
    elseif nargin == 3
        ax = NaN;
        raster = varargin{1};
        histogram = varargin{2};
        edges = varargin{3};
    elseif nargin == 1 || nargin == 2;
        if nargin == 2
            ax = varargin{1};
            raster = varargin{2};
        else
            ax = NaN;
            raster = varargin{1};
        end
        
        bigraster = [];
        
        for ii = 1:numel(raster)
            line = raster{ii};
            bigraster = [bigraster; reshape(line,numel(line),1)];
        end
        
        edges = linspace(min(bigraster),max(bigraster),100);
        histogram = histc(bigraster,edges);
    else
        error 'What the shit is this gimme something to work with c''mon';
    end
    
    if isnan(ax)
        ax = gca;
    end
    
    hold on;
    bar(ax,edges,histogram,'histc');
    
    maxCount = max(histogram);
    nLines = numel(raster);
    
    if maxCount == 0
        lineHeights = 1:nLines;
    else
        lineHeights = maxCount + (1:nLines)*maxCount/nLines;
    end
    
    for ii = 1:nLines;
        line = raster{ii};
        height = lineHeights(ii);
        
        plot(ax,line,height*ones(size(line)),'LineStyle','none','Marker','.','Color','k');
    end
    
    xlim([edges(1) edges(end)]);
    
    if maxCount > 0
        ylim([0 2*maxCount]);
    end
end