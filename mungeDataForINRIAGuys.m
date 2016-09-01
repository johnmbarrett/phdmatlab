function mungeDataForINRIAGuys(recording)
    [fileDir,filename] = getAnalysisOutputDir(recording);
    load([fileDir '\' filename '_photodiode_timings.mat']);
    load([fileDir '\' filename '_psrh_abs.mat']);
    load([fileDir '\' filename '_spikes_concat_forceclustered_no_ignorenoise_no.mat']);
    
    %%
    channels = str2double(cellstr(char(cells(:,1:2))));
    clusters = cells(:,3);

    %%
    widths = valuess{1};
    phases = valuess{2};
    [sortedWidths,widthIndices] = sort(widths);
    [sortedPhases,phaseIndices] = sort(phases);
    sortedRasters = squeeze(rasters(:,widthIndices,phaseIndices,1));
    sortedRepeats = squeeze(repeats(:,widthIndices,phaseIndices,1));
    sortedTimings = squeeze(stimulusTimings(widthIndices,phaseIndices,1));
    sortedRepeats = squeeze(sortedRepeats(1,:,:));
    N = sum(sum(sortedRepeats));
    allSpikes = cell(size(rasters,1)*N,1);
    %%
    nn = 0;

    for ii = 1:size(rasters,1)
        for jj = 1:numel(widths)
            for kk = 1:numel(phases)
                for ll = 1:sortedRepeats(jj,kk)
                    nn = nn + 1;
                    allSpikesPrefix = [ii numel(phases)*(jj-1)+kk ll];
                    raster = sortedRasters{ii,jj,kk};
                    line = raster{ll};

                    if isempty(line)
                        allSpikes{nn} = allSpikesPrefix;
                        continue;
                    end

                    presentationTime = sortedTimings{jj,kk}(2,2,ll)-sortedTimings{jj,kk}(1,1,ll);
                    line = line(line >= 0 & line <= presentationTime);

                    if isempty(line)
                        allSpikes{nn} = allSpikesPrefix;
                        continue;
                    end

                    allSpikes{nn} = [allSpikesPrefix reshape(line,1,numel(line))];
                end
            end
        end
    end
    %%
    stimuli = zeros(numel(widths)*numel(phases),2);
    for jj = 1:numel(widths)
        for kk = 1:numel(phases)
            stimuli(numel(phases)*(jj-1)+kk,:) = [25*sortedWidths(jj)/2 8*sortedPhases(kk)];
        end
    end
    %%
    save([filename '_all_spikes.mat'],'allSpikes','channels','clusters','stimuli');
end