function nexcat(outputFile)
    nexFiles = uipickfiles('FilterSpec','*.nex');
    
    assert(numel(nexFiles) > 0);
    
    if numel(nexFiles) == 1
        disp('Only one nex file selected, nothing more to do.');
        return;
    end
    
    allNames = {};
    allTimestamps = {};
    allWaveforms = {};
    msgFn = @(x) ['Processing nex file #' num2str(x) ' of ' num2str(numel(nexFiles))];
    wf = waitbar(0,msgFn(1));
    
    maxT = 0;
    for ii = 1:numel(nexFiles)
        if ii > 1
            waitbar((ii-1)/numel(nexFiles),wf,msgFn(ii));
        end
        
        nex = readNexFile(nexFiles{ii});
        
        channelFn = @(x) ['Processing channel #' num2str(x) ' of ' num2str(numel(nex.waves))];
        wc = waitbar(0,channelFn(1));
        
        for jj = 1:numel(nex.waves)
            if jj > 1
                waitbar((jj-1)/numel(nex.waves),wc,channelFn(jj));
            end
            
            wave = nex.waves{jj};
            name = wave.name;
            index = find(strcmp(name,allNames));
            
            if isempty(index)
                allNames{end+1} = name; %#ok<*AGROW>
                allTimestamps{end+1} = double(wave.timestamps)+maxT;
                allWaveforms{end+1} = double(wave.waveforms);
            else
                allTimestamps{index} = [allTimestamps{index}; double(wave.timestamps)+maxT];
                allWaveforms{index} = [allWaveforms{index} double(wave.waveforms)];
            end
        end
        
        maxT = maxT + nex.tend;
        
        close(wc);
    end
    close(wf);
            
    bigNex = nexCreateFileData(nex.freq);
    
    ws = waitbar(0,'Saving Nex file...');
    for ii = 1:numel(allNames)
        bigNex = nexAddWaveform(bigNex,bigNex.freq,allTimestamps{ii},allWaveforms{ii},allNames{ii});
        waitbar(ii/numel(allNames),ws);
    end
    
    close(ws);
    
%     maxT = -Inf;
%     
%     msgFn = @(x) ['Processing nex file #' num2str(x) ' of ' num2str(numel(nexFiles))];
%     
%     
%     for ii = 1:numel(nex1.waves)
%         wave = nex1.waves{ii};
%         bigNex = nexAddWaveform(bigNex,nex1.freq,wave.timestamps,wave.waveforms,wave.name);
%         maxT = max(maxT,max(wave.timestamps));
%     end
%     
%     wf = waitbar(1/numel(nexFiles),wf,msgFn(1));
%     
%     for ii = 2:numel(nexFiles);
%         oldMaxT = maxT;
%         nex = readNexFile(nexFiles{ii});
%     
%         channelFn = @(x) ['Processing channel #' num2str(x) ' of ' num2str(numel(nex.waves))];
%         wc = waitbar(0,channelFn(1));
%         
%         for jj = 1:numel(nex.waves);
%             if jj > 1
%                 waitbar((jj-1)/numel(nex.waves),wc,channelFn(jj));
%             end
%             
%             waves = vertcat(bigNex.waves{:});
%             names = {waves.name};
%             wave = nex.waves{ii};
%             name = wave.name;
%             index = find(strcmp(name,names));
%             
%             if isempty(index)
%                 bigNex = nexAddWaveform(bigNex,nex.freq,double(wave.timestamps+oldMaxT),double(wave.waveforms),wave.name);
%             else
%                 bigNex.waves{index}.timestamps = [bigNex.waves{index}.timestamps; double(wave.timestamps+oldMaxT)];
%                 bigNex.waves{index}.waveforms = [bigNex.waves{index}.waveforms double(wave.waveforms)];
%             end
%             
%             maxT = max(maxT,max(wave.timestamps+oldMaxT));
%         end
%         
%         close(wc);
%         wf = waitbar(ii/numel(nexFiles),wf,msgFn(ii));
%     end
%     
%     close(wf);
    
    writeNexFile(bigNex,outputFile);
end