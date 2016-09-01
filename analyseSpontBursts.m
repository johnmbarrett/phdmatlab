function analyseSpontBursts(recording,varargin)
    opt = getopt('bursts=''no'' frs=''no'' overwrite=''no'' window=60',varargin{:});
    window = opt.window;
    
    [allSpikeTimes,channels] = getAllSpikeTimes(recording,true,true);
    maxT = ceil(max(vertcat(allSpikeTimes{:}))/window);
    nCells = numel(allSpikeTimes);
    cellsPerWindow = zeros(maxT,1);
    windowedSpikeTimes = cell(maxT,nCells);
    firingRates = zeros(maxT,nCells);
    
    for ii = 1:nCells
        spikeTimes = allSpikeTimes{ii};
        
        for jj = 1:maxT
            isInWindow = spikeTimes >= (jj-1)*window & spikeTimes < jj*window;
            
            if sum(isInWindow) > 1
                cellsPerWindow(jj) = cellsPerWindow(jj)+1;
            end
            
            windowedSpikeTimes{jj,ii} = spikeTimes(isInWindow);
            firingRates(jj,ii) = sum(isInWindow);
        end
    end
    
    labels = num2str(channelIndexToMCSChannelNumber(1:60)');
    
    function channelIndices = getChannelIndices(channelIndex,single)
        if single
            set(gca,'XTick',[],'YTick',[]);
            label = labels(channelIndex,:);
            channelIndices = strcmp(label,cellstr(channels));
        else
            channelIndices = 1:size(firingRates,2);
        end
    end

    x = (1:maxT)*window/60;
    
    function plotFn(datas,channelIndex,~,~)
        channelIndices = getChannelIndices(channelIndex,nargin < 4);

        if isempty(channelIndices) || sum(channelIndices) == 0
            return;
        end
        
        if iscell(datas)
            prctiles = zeros(maxT,3);
        
            for kk = 1:maxT
                data = vertcat(datas{kk,channelIndices});
                
                if isempty(data)
                    continue;
                end
                
                prctiles(kk,:) = prctile(data,[25 50 75]);
            end
        elseif isnumeric(datas)
            if sum(channelIndices) == 1
                plot(x,datas(:,channelIndices));
                
                if nargin < 4
                    set(gca,'XTick',[],'YTick',[]);
                end
                
                return;
            end

            prctiles = prctile(datas(:,channelIndices),[25 50 75],2);
        end
        
        y = prctiles(:,2);
        e = abs(prctiles(:,[1 3])-repmat(y,1,2));

        boundedline(x,y,e);
    end
    
    plotFnFn = @(data) @(channelIndex,varargin) plotFn(data,channelIndex,varargin{:});
    
    fileDir = getAnalysisOutputDir(recording);
    
    function doPlots(plotFn,meaSuffix,allSuffix,ylab)
        if ~isnan(meaSuffix)
            mcRackPlot(plotFn);
            subplot(8,8,7*8+1);
            set(gcf,'Position',[0 0 1600 900]);
            set(gca,'XTick',[],'YTick',[]);
            xlabel('Time');
            ylabel(ylab);

            figfile = sprintf('%s\\%s_%s',fileDir,fileDir,meaSuffix);
            saveas(gcf,figfile,'fig');
            saveas(gcf,figfile,'png');
        end
        
        figure;
        plotFn(NaN,NaN,NaN); % the best code
        xlabel('Time/min');
        ylabel(ylab);
        title(recording);

        figfile = sprintf('%s\\%s_%s',fileDir,fileDir,allSuffix);
        saveas(gcf,figfile,'fig');
        saveas(gcf,figfile,'png');
    end

    doPlots(@(varargin) plot(x,cellsPerWindow),NaN,'ncells','# Cells Firing');
    
    if strcmpi(opt.frs,'yes')
        firingRates = firingRates/window;

        doPlots(plotFnFn(firingRates),'fr_mea','fr_all','Firing Rate/Hz');
    end

    if strcmpi(opt.bursts,'yes')
        burstfile = sprintf('%s\\%s_burst_analysis.mat',fileDir,fileDir);
        
        if strcmpi(opt.overwrite,'yes') || ~exist(burstfile,'file')
            burstiness = zeros(maxT,nCells); 
            burstLens = cell(maxT,nCells);
            burstDurs = cell(maxT,nCells);
            ibis = cell(maxT,nCells);

            for ii = 1:maxT
               for jj = 1:nCells
                   t = windowedSpikeTimes{ii,jj};

                   if numel(t) < 2
                       continue;
                   end

                   [~,len,starts] = burstold(t);

                   if isempty(len)
                       continue;
                   end

                   burstiness(ii,jj) = sum(len)/numel(t);
                   burstLens{ii,jj} = len(:);
                   burstDurs{ii,jj} = diff([t(starts) t(starts+len-1)],[],2);
                   ibis{ii,jj} = diff(t(starts));
               end
            end

            save(burstfile,'burstiness','burstLens','burstDurs','ibis');
        else
            load(burstfile);
        end
       
       datasets = {burstiness burstLens burstDurs ibis};
       suffixes = {'burstiness' 'burst_len' 'burst_dur' 'ibi'};
       ylabs = {'Burstiness (a.u.)' 'Burst Length (# spikes)' 'Burst Duration/s' 'Inter-Burst Interval/s'};
       
       for ii = 1:numel(datasets)
           doPlots(plotFnFn(datasets{ii}),sprintf('%s_mea',suffixes{ii}),sprintf('%s_all',suffixes{ii}),ylabs{ii});
       end
    end
end