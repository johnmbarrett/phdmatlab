function [nTriggers,nTriggerss,stimSpecs,specIndices,stimFileIndices] = parseULEDArrayBatchFiles(filePaths)
    if ~iscell(filePaths)
        filePaths = {filePaths};
    end
    
    nBatches = numel(filePaths);
    stimFiles = {};
    stimFileIndices = cell(nBatches,1);
    
    for ii = 1:nBatches
        filePath = filePaths{ii};
        [~,fileName,extension] = fileparts(filePath);
    
        if strcmp(extension(2:end),'txt')
            stimFiles{end+1} = filePath; %#ok<AGROW>
        elseif strcmp(extension(2:end),'bat')
            fin = fopen(filePath,'r');
            
            fscanf(fin,'uledctrl.exe COM%d %d');
            
            while ~feof(fin)
                stimFile = textscan(fin,' %q');
                nFiles = numel(stimFile{1});
                
                if ii == 1
                    stimFileIndices{ii} = (1:nFiles);
                else
                    stimFileIndices{ii} = stimFileIndices{ii-1}(end)+(1:nFiles);
                end
                
                stimFiles(end+(1:nFiles)) = stimFile{1}; %#ok<AGROW>
            end
        else
            warning('Unrecognised format for file %d (%s), skipping...\n',ii,fileName);
        end
    end
    
    nStimFiles = numel(stimFiles);
    nTriggerss = zeros(nStimFiles,1);
    specIndices = zeros(nStimFiles,1);
    stimSpecs = {};
    
    for ii = 1:nStimFiles
        stimFile = stimFiles{ii};
        
        seen = find(strcmp(stimFile,stimFiles),1);
        
        if ~isempty(seen) && seen < ii;
            nTriggerss(ii) = nTriggerss(seen);
            specIndices(ii) = specIndices(seen);
            continue;
        end
        
        fin = fopen(stimFile,'r');
        nTriggers = 0;
        stimSpec = struct('index',{},'type',{},'X',{},'Y',{},'width',{},'height',{},'period',{},'duration',{},'minRadius',{},'maxRadius',{},'path',{},'direction',{},'squareSize',{},'seed',{},'p',{});
        n = 0;
        
        while ~feof(fin)
            s = fgetl(fin);
            n = n + 1;
            stimSpec(n).type = s(1);
            stimSpec(n).index = nTriggers + 1;
            
            switch s(1)
                case 'b'
                    A = sscanf(s,'b %d');
                    stimSpec(n).duration = A;
                    nTriggers = nTriggers + 1;                    
                case 'c'
                    A = sscanf(s,'c %d %d %d %d %d %d');
                    stimSpec(n).type = e;
                    stimSpec(n).direction = 'in';
                    stimSpec(n).X = A(1);
                    stimSpec(n).Y = A(2);
                    stimSpec(n).minRadius = A(3);
                    stimSpec(n).maxRadius = A(4);
                    stimSpec(n).width = A(5);
                    stimSpec(n).period = A(6);
                    nTriggers = nTriggers + stimSpec(n).maxRadius - stimSpec(n).minRadius + 1;
                case 'e'
                    A = sscanf(s,'c %d %d %d %d %d %d');
                    stimSpec(n).type = e;
                    stimSpec(n).direction = 'out';
                    stimSpec(n).X = A(1);
                    stimSpec(n).Y = A(2);
                    stimSpec(n).minRadius = A(3);
                    stimSpec(n).maxRadius = A(4);
                    stimSpec(n).width = A(5);
                    stimSpec(n).period = A(6);
                    nTriggers = nTriggers + stimSpec(n).maxRadius - stimSpec(n).minRadius + 1;
                case 'g'
                    gifPath = s(3:end);
                    stimSpec(n).path = gifPath;
                    
                    I = imread(gifPath,'gif');
                    nFrames = size(I,4);
                    nTriggers = nTriggers + nFrames;
                    
                    gifInfo = imfinfo(gifPath,'gif');
                    delayTime = gifInfo.DelayTime*10;
                    
                    stimSpec(n).period = delayTime;
                    stimSpec(n).duration = delayTime*nFrames;
                case 'i'
                    A = sscanf(s,'i %d %s');
                    stimSpec(n).duration = A(1);
                    stimSpec(n).path = char(A(2:end))';
                    nTriggers = nTriggers + 1;                    
                case 'm'
                    A = sscanf(s,'m %s %d %d');
                    stimSpec(n).direction = char(A(1:end-2))';
                    stimSpec(n).period = A(end-1);
                    stimSpec(n).width = A(end);
                    nTriggers = nTriggers + 16 - stimSpec(n).width;
                case 'r'
                    A = sscanf(s,'r %d %d %d %d %d');
                    stimSpec(n).X = A(1);
                    stimSpec(n).Y = A(2);
                    stimSpec(n).width = A(3);
                    stimSpec(n).height = A(4);
                    stimSpec(n).duration = A(5);
                    nTriggers = nTriggers + 1;                    
                case 'w'
                    A = sscanf(s,'w %d %d %d %d %d %d %d %d %d');
                    stimSpec(n).X = A(1);
                    stimSpec(n).Y = A(2);
                    stimSpec(n).width = A(3);
                    stimSpec(n).height = A(4);
                    stimSpec(n).squareSize = A(5);
                    stimSpec(n).period = A(6);
                    stimSpec(n).duration = A(7)*stimSpec(n).period;
                    stimSpec(n).seed = A(8);
                    stimSpec(n).p = A(9);
                    nTriggers = nTriggers + A(7);
            end
        end
        
        nTriggerss(ii) = nTriggers;
        stimSpecs{end+1} = stimSpec; %#ok<AGROW>
        specIndices(ii) = numel(stimSpecs);
    end
    
    nTriggers = sum(nTriggerss);
end
