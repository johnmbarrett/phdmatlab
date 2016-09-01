function allCellsRaster(recording,time,varargin)
    options = getopt('axes filesuffix=''_MCD_trimmed_spikes'' isclustered=1',varargin{:});
    fileSuffix = options.filesuffix;
    
    assert(isscalar(options.isclustered) && ismember(options.isclustered,[-1 0 1]),'isclustered must be one of -1 (force unclustered), 0 (clustered if availabel), 1 (force clustered)');
    isClustered = options.isclustered;
    
    if isempty(options.axes) || isnan(options.axes)
        figure;
        ax = axes;
    else
        ax = options.axes;
    end
    
    hold on;
    
    [fileDir,filename] = getAnalysisOutputDir(recording);
    
    if nargin < 2 || isempty(time) || any(~isfinite(time(:)))
        [err,file] = ns_OpenFile([filename '.mcd']);
        
        if err
            time = [0 Inf];
        else
            [err,fileInfo] = ns_GetFileInfo(file);
            
            if err
                time = [0 Inf];
            else
                time = [0 fileInfo.TimeSpan];
            end
            
            ns_CloseFile(file);
        end
    elseif isscalar(time)
        time = [0 time];
    end
    
    if nargin < 1
        error 'What, precisely, do you expect me to do here?';
    end
    
    row = 1;
    ticks = [];
    labels = {};
    
    files = dir([fileDir '\*' filename '_channel_*' fileSuffix '.mat']);
    filenames = {files.name};
    rawFiles = cellfun('isempty',strfind(filenames,'\times_'));
    
    if isClustered == 1 || (isClustered == 0 && ~all(rawFiles))
        filenames = filenames(~rawFiles);
        usingClustered = true;
    else
        filenames = filenames(rawFiles);
        usingClustered = false;
    end
    
    nFiles = numel(filenames);
    
    for channel = 1:nFiles
        file = [fileDir '\' filenames{channel}];
        
        if usingClustered
            load(file,'cluster_class');
        else
            load(file,'index');
            cluster_class = [ones(numel(index),1) index'];
        end
        
        if isempty(cluster_class(:,2));
            continue;
        end
        
        clusters = cluster_class(:,1);
        
        nClusters = max(clusters);
        
        if nClusters == 0
            continue;
        end
        
        if isempty(ticks)
            ticks = nClusters;
        else
            ticks(end+1) = ticks(end) + nClusters; %#ok<AGROW>
        end
        
        label = regexp(file,['channel_(.+)' fileSuffix],'tokens');
        labels(end+1) = label{1}; %#ok<AGROW>
        
        spikeTimes = cluster_class(:,2)/100;
        
        for cluster = 1:nClusters
            spikes = spikeTimes(spikeTimes > time(1) & spikeTimes <= time(2) & clusters == cluster)';
            
            if isempty(spikes)
                continue;
            end
            
            line(repmat(spikes,2,1),[row;row+1]*ones(size(spikes)),'Color','k'); %,'LineStyle','none','Marker','.');
            row = row+1;
        end
    end
    
    set(ax,'YTick',ticks);
    set(ax,'YTickLabel',labels);
    
    xlim(ax,time);
    xlabel('Time/s');
    ylim(ax,[0 row+1]);
    ylabel('Channel');
end
            