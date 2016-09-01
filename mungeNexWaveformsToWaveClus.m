function mungeNexWaveformsToWaveClus(recording,varargin)
    options = getopt('nexfiles={} choosenexfiles=false sorted=true prefixlength=4 filesuffix=''_MCD_trimmed_spikes.mat''',varargin{:});
    
    [fileDir,filename] = getAnalysisOutputDir(recording);

    if all(logical(options.choosenexfiles))
        nexFiles = uipickfiles('FilterSpec','*.nex');
    elseif isempty(options.nexfiles)
        nexFiles = {[filename ' spikes trimmed sorted.nex']};
    elseif ischar(options.nexfiles)
        nexFiles = {options.nexfiles};
    elseif iscell(options.nexfiles)
        nexFiles = options.nexfiles;
    else
        warning('No valid .nex files specified for recording %s, nothing to do.',filename);
        return;
    end
    
    isSorted = all(logical(options.sorted));
        
    for hh = 1:numel(nexFiles)
        nexFile = nexFiles{hh};
        
        if isempty(filename)
            [~,outFile] = fileparts(nexFile);
            outDir = outFile;
        else
            outFile = filename;
            outDir = fileDir;
        end

        if ~ischar(nexFile)
            warning('File #%d for recording %s is not a valid file path, ignoring...',hh,filename);
            continue;
        elseif ~exist(nexFile,'file') || ~strcmp(nexFile(end-3:end),'.nex')
            warning('File path %s for recording %s does not point to a valid .nex file, ignoring...',nexFile,filename);
            continue;
        end
        
        nex = readNexFile(nexFile);
        waves = vertcat(nex.waves{:});
        waveNames = {waves.name}';
        
        if isSorted
            % need this step if exported templates as single-unit waveforms
            waveNames = vertcat(waveNames{~cellfun('isempty',strfind(waveNames,'wf'))});
        else
            waveNames = vertcat(waveNames{:});
        end
        
        channelLabels = unique(waveNames(:,1:options.prefixlength),'rows');

        for ii = 1:size(channelLabels,1)
            channelLabel = channelLabels(ii,:);

            units = waves(strncmp(channelLabel,cellstr(waveNames),options.prefixlength));
            allWaves = [units.waveforms]';

            spikeTimes = vertcat(units.timestamps)*100;
           
            if ~isSorted
                index = spikeTimes'; %#ok<NASGU>
                spikes = allWaves; %#ok<NASGU>
                save([outDir '\' outFile '_channel_' channelLabel(3:options.prefixlength) options.filesuffix],'spikes','index');
                continue;
            end
                
            cluster_class(:,2) = spikeTimes;

            wavesSoFar = 0;
            for jj = 1:numel(units)
                unitLetter = units(jj).name(options.prefixlength+1);

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

            save([outDir '\times_' outFile '_channel_' channelLabel(3:options.prefixlength) '_MCD_trimmed_spikes.mat'],'spikes','cluster_class');
            clear cluster_class;
        end
    end
end