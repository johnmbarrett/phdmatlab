function mcRackPlot(data,file,x,yy)
    if nargin < 1
        error('What, exactly, would you like me to plot?');
    end
    
    isPlotFn = isa(data,'function_handle');
    
    if ~isPlotFn
        if nargin < 4
            yy = [-100 100];
        end

        if nargin < 3
            x = (1:size(data,1))/file.sampleRate;
        end
    end
    
    if nargin < 2 || isempty(file)
        file = struct('electrodes',struct('label',cellstr(num2str(channelIndexToMCSChannelNumber(1:60)'))));
    end
    
    % TODO : refactor
    figure;
    
    for ii = 1:60
        chlabel = str2double(file.electrodes(ii).label);
        subplot(8,8,getMCSSubplotIndex(chlabel));
        
        if isPlotFn
            data(ii,gca);
            continue;
        end
        
        plot(x,data(:,ii),'Color',[0 1 0]);
        set(gca,'Color',[0 0 0],'XTick',[],'YTick',[]);
        
        if size(yy,1) == 1
            ylim(yy);
        else
            ylim(yy(ii,:));
        end
    end
end