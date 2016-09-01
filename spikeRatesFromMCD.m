%%
ns_SetLibrary('D:\John B\nsMCDlibrary.dll');
filenames = {'Spont 1', 'Stim 1', 'Spont 2', 'Stim 2', 'Spont 3', 'Stim 3', 'Spont 4'};
rates = zeros(60,7);
spikeTimess = cell(60,7);
recordingTimes = zeros(7,1);
channelLabels = zeros(60,1);

for ii = 1:7
    filename = filenames{ii};
    
    [error,file] = ns_OpenFile([filename '.mcd']);
    
    if error
        fprintf('Could not open file %s.mcd (%d/7)\n',filename,ii);
        continue;
    end
    
    [error,fileInfo] = ns_GetFileInfo(file);
    
    if error
        fprintf('Could not load file info from file %s.mcd (%d/7)\n',filename,ii);
        continue;
    end
    
    [error,entities] = ns_GetEntityInfo(file,1:fileInfo.EntityCount);
    
    if error
        fprintf('Could not load entity info from file %s.mcd (%d/7)\n',filename,ii);
        continue;
    end
    
    recordingTimes(ii) = fileInfo.TimeSpan;
    
    for jj = 1:60
        entity = entities(jj);
        channelLabel = str2double(entity.EntityLabel(end-1:end));
        
        if ii == 1
            kk = jj;
            channelLabels(jj) = channelLabel;
        else
            kk = find(channelLabels == channelLabel);
        end
        
        if numel(kk) ~= 1
            fprintf('Unknown channel label for channel %d of file %d, kk = \n',jj,ii);
            disp(kk);
            continue;
        end
        
        rates(kk,ii) = entity.ItemCount/fileInfo.TimeSpan;
        
        if entity.ItemCount > 0
            [error,spikeTimes] = ns_GetSegmentData(file,jj,1:entity.ItemCount);
            
            if error
                fprintf('Could not load spike times for channel %d of file %d',jj,ii);
                continue;
            end
            
            spikeTimess{kk,ii} = spikeTimes;
        end
    end
end

%%
figure
for ii = 1:60
    channelLabel = channelLabels(ii);
    row = mod(channelLabel,10)-1;
    col = floor(channelLabel/10);
    plotnum = 8*row+col;
%     fprintf('%02d %02d %02d %02d %02d\n',ii,channelLabel,row,col,plotnum);
    subplot(8,8,plotnum);
    bar(1:7,rates(ii,:));
end

%%
colours = 'rgbmyck';
figure
for ii = 1:60
    channelLabel = channelLabels(ii);
    row = mod(channelLabel,10)-1;
    col = floor(channelLabel/10);
    plotnum = 8*row+col;
    subplot(8,8,plotnum);
    hold on;
    
    maxRate = -Inf;
    
    for jj = 1:7
        spikeTimes = spikeTimess{ii,jj};
        
        if isempty(spikeTimes)
            continue;
        end
        
        edges = 0:ceil(recordingTimes(jj));
        h = histc(spikeTimes,edges);
        maxRate = max(maxRate,max(h));
        plot(edges,h,colours(jj));
    end
    
    if maxRate > 0
        ylim([0 maxRate]);
    end
    
    xlim([0 max(recordingTimes)]);
end