function mungeNexWaveformsToWaveClusForHarbi(nexFiles,outputDir,clustered)
    if nargin < 3
        clustered = false;
    end
    
    if nargin < 2
        outputDir = '.';
    end
    
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end

    if nargin < 1
        nexFiles = uipickfiles('FilterSpec','*.nex');
    elseif ischar(nexFiles)
        nexFiles = {nexFiles};
    elseif ~iscell(nexFiles)
        error('Input must be a cell array of paths to .nex files');
    end
        
    for hh = 1:numel(nexFiles)
        nexFile = nexFiles{hh};
        [~,nexFilename] = fileparts(nexFile);
        
        if ~exist(sprintf('%s/%s',outputDir,nexFilename),'dir')
            mkdir(sprintf('%s/%s',outputDir,nexFilename));
        end

        if ~ischar(nexFile)
            warning('File #%d for recording %s is not a valid file path, ignoring...',hh,recording.dataFile);
            continue;
        elseif ~exist(nexFile,'file') || ~strcmp(nexFile(end-3:end),'.nex')
            warning('File path %s for recording %s does not point to a valid .nex file, ignoring...',nexFile,recording.dataFile);
            continue;
        end
        
        nex = readNexFile(nexFile);
        waves = vertcat(nex.waves{:});
        waveNames = vertcat(waves.name);
        channelLabels = unique(waveNames(:,1:4),'rows');

        for ii = 1:size(channelLabels,1)
            channelLabel = channelLabels(ii,:);

            units = waves(strncmp(channelLabel,cellstr(waveNames),4));
            allWaves = [units.waveforms]';

            if clustered
                cluster_class(:,2) = vertcat(units.timestamps)*100;

                wavesSoFar = 0;
                for jj = 1:numel(units)
                    unitLetter = units(jj).name(5);

                    if unitLetter == 'U'
                        cluster = 0;
                    else
                        cluster = (unitLetter - 'a') + 1;
                    end

                    nWaves = numel(units(jj).timestamps);
                    cluster_class(wavesSoFar+(1:nWaves),1) = cluster; %#ok<AGROW>
                    wavesSoFar = wavesSoFar+nWaves;
                end

                [cluster_class,sortIndices] = sortrows(cluster_class,2);
                spikes = allWaves(sortIndices,:); %#ok<NASGU>
                save(sprintf('%s/%s/%s.mat',outputDir,nexFilename,channelLabel),'spikes','cluster_class');
            else
                spikes = allWaves;
                index = vercat(units.timestampes)*100;
                
                save(sprintf('%s/%s/%s.mat',outputDir,nexFilename,channelLabel),'spikes','index');
            end
            
            clear cluster_class;
        end
    end
end