function mungeCStimOutputIntoStimulateOutput(stimulusFile,outputFile)
    fin = fopen(stimulusFile,'r');
    closeFin = onCleanup(@() fclose(fin));
    
    reg = '^cstimx?p?.exe -r ([0-9]+) -g ([0-9]+) -b ([0-9]+) -p ([0-9]+) -t ([0-9]+) -o "?([^"]+)"?$';
    
    getExtraParams = {};
    stimuli = [];
    vbls = [];
    ii = 0;
    
    while ~feof(fin)
        [~,~,~,~,tokens] = regexp(fgetl(fin),reg);
        
        if isempty(tokens)
            continue;
        end
        
        ii = ii + 1;
        
        tokens = tokens{1};
        timingFilePath = tokens{6};
        tokens = str2double(tokens(1:5));
        
        colour = tokens(1:3);
        period = floor(60*tokens(4)/1000);
        maxT = tokens(5)*60;
        
        getExtraParams{ii} = struct('colour',colour,'period',period,'maxT',maxT); %#ok<AGROW>
        
        slashes = strfind(timingFilePath,'\');
        currentPath = pwd;
        
        if isempty(slashes)
            timingPath = pwd;
            timingFile = timingFilePath;
        else
            lastSlash = slashes(end);
            timingPath = timingFilePath(1:lastSlash);
            timingFile = timingFilePath(lastSlash+1:end);
        end
        
        try
            cd(timingPath)
        catch err
            warning('Error "%s" while trying to CD to timing file directory, looking in current directory for a copy...',err.message); %#ok<WNTAG>
        end
        
        try
            tin = fopen(timingFile,'r');
        catch err
            warning('Error "%s" trying to open timing file %s, skipping...',err.message,timingFile); %#ok<WNTAG>
            continue;
        end
        
        if tin < 0
            warning('Failed to open timing file %s, skipping...',timingFile); %#ok<WNTAG>
            continue;
        end
        
        closeTin = onCleanup(@() fclose(tin));
        
        % TODO : option to choose timing source
        systemStartTime = sscanf(fgetl(tin),'%d\t%d\t%d\t%d\t%d\t%d\t%d\n');
        systemStartTime = datenum([systemStartTime(1:5); systemStartTime(6)+systemStartTime(7)/1000]')*24*3600;
        
        fileStartTime = sscanf(fgetl(tin),'File Time:\t%d\t%d\n');
        fileStartTime = (fileStartTime(1)*10^(numel(num2str(fileStartTime(2))))+fileStartTime(2))/1e7; %#ok<NASGU>
        
        tickStartTime = sscanf(fgetl(tin),'Tick Count:\t%d\n');
        tickStartTime = tickStartTime/1000;
        
        unbiasedStartTime = sscanf(fgetl(tin),'Unbiased:\t%d\n');
        unbiasedStartTime = unbiasedStartTime/1e7; %#ok<NASGU>
        
        performanceFrequency = sscanf(fgetl(tin),'Perf Freq:\t%d\n');
        
        performanceStartTime = sscanf(fgetl(tin),'Perf Count:\t%ld\n');
        performanceStartTime = performanceStartTime/performanceFrequency; %#ok<NASGU>
        
        if ii == 1
            timeOffset = systemStartTime-tickStartTime; %#ok<NASGU>
        end
        
        sscanf(fgetl(tin),'---\n');
        
        while ~feof(tin)
            ts = sscanf(fgetl(tin),'%d\t%d\t%d\t%d\n');
            stimuli = [stimuli; ii]; %#ok<AGROW>
            vbls = [vbls; ts(2)/1e3]; %#ok<AGROW>
        end
        
        cd(currentPath);
    end
    
    save(outputFile,'getExtraParams','stimuli','timeOffset','vbls');
end