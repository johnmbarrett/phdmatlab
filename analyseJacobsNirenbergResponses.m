function analyseJacobsNirenbergResponses(recording,varargin)
    %% INTITIALIZATION and defining significant responses
    bp = defaultParams;
    bp.samp_iter = 100;
    
    fileDir = getAnalysisOutputDir(recording);
    
    opts = getopt('latencytype=''firstspike'' significance=''rf'' usesavedposteriors=false ifrfit=''gauss'' ifrparam=1 latencyparam=50 posteriorsonly=false seed=NaN',varargin{:});
    
    if isstruct(opts.seed)
        rng(opts.seed);
    end
    
    if ~opts.usesavedposteriors
        load([fileDir '\' recording.dataFile '_psrh_abs.mat']);

        try
            load([fileDir '\' recording.dataFile '_spikes_concat_forceclustered_yes_ignorenoise_no.mat']);
        catch %#ok<CTCH>
            load([fileDir '\' recording.dataFile '_spikes_concat_forceclustered_no_ignorenoise_no.mat']);
        end

        load(recording.dataFile,'getExtraParams');
        
        % TODO : use this later to retreive stimulus params
        stimulusParamss = vertcat(getExtraParams{2:2:end}); %#ok<USENS>
        stimulusParamss = [stimulusParamss.(factors{1}); stimulusParamss.(factors{2})]'.*repmat([25/2 8],size(stimulusParamss,1),1); %#ok<USENS>
        [~,sortedStimulusIndices] = sortrows(stimulusParamss);
        
        load([fileDir '\' recording.dataFile '_photodiode_timings']);

        stimulusTimes = stimulusTimes - recordingStartTime; %#ok<NODEF>

        nStimuli = floor((numel(stimulusTimes)-1)/2);
        gratings = 2*(1:nStimuli);
        masks = gratings+1;
        maxStimulusLengths(1) = max(stimulusTimes(masks)-stimulusTimes(gratings));
        maxStimulusLengths(2) = max(stimulusTimes(masks+1)-stimulusTimes(masks));

        widths = valuess{1}*25/2; %#ok<USENS>
        phases = valuess{2}*8;

    %     responsiveChannels = [13 2; 17 1; 35 2; 36 1; 37 1; 37 2; 48 2; 58 1; 62 1; 71 2; 73 1];

        if strcmp(opts.significance,'rf')
            recordings = initRecordings;
            if isfield(recording,'rfFileIndex')
                rfRecording = recordings(recording.rfFileIndex);

                significantResponses = false(nCells,1);

                for ii = 1:nCells
                    rfDir = getAnalysisOutputDir(rfRecording);
                    rfFile = [rfDir '\srrf_' rfRecording.dataFile '_channel_' char(cells(ii,1:2)) '_cluster_' num2str(cells(ii,3)) '_clustered_phototrigger_all_spikes.mat']; %#ok<FDEPR>

                    if exist(rfFile,'file')
                        load(rfFile);
                        significantResponses(ii) = logical(significant);
                    end
                end
            else
                significantResponses = true(nCells,1);
            end

            responsiveCells = find(significantResponses);
        elseif strcmp(opts.significance,'visual')
            lightResponses = importdata([fileDir '\lightresponses.txt'],'\t');
            lightResponses = [double(num2str(lightResponses(:,1))) lightResponses(:,2)+1];
            responsiveCells = cellIndices(sub2ind(size(cellIndices),lightResponses(:,1),lightResponses(:,2),lightResponses(:,3)));
        elseif strcmp(opts.significance,'debug')
            lightResponses = [double(num2str(41)) 1];
            responsiveCells = cellIndices(sub2ind(size(cellIndices),lightResponses(:,1),lightResponses(:,2),lightResponses(:,3)));
        else
            responsiveCells = find(cells(:,3) > 0); %#ok<NODEF>
        end

        nResponses = numel(responsiveCells);

        cellSpikeTimes = cell(nResponses,1);

        for ii = 1:nResponses
            cellIndex = responsiveCells(ii);

            channel = cells(cellIndex,1:2);
            cluster = cells(cellIndex,3); 

            cellSpikeTimes{ii} = sortedSpikeTimes( ...
                sortedSpikeChannels(:,1) == channel(1) & ...
                sortedSpikeChannels(:,2) == channel(2) & ...
                sortedSpikeClusters == cluster); %#ok<NODEF>
        end

        [widths,widthIndices] = sort(widths); %#ok<NASGU>
        [phases,phaseIndices] = sort(phases); %#ok<NASGU>
        
%         if false

        %     sortedRasters = squeeze(rasters(responsiveCells,widthIndices,shiftIndices,1));
        %     sortedTimings = squeeze(stimulusTimings(widthIndices,shiftIndices,1));
            binWidth = 1/1000;
            binss = {(0:ceil(maxStimulusLengths(1)/binWidth))' (0:ceil(maxStimulusLengths(2)/binWidth))'};
            psths = {zeros(numel(binss{1}),8,4,nResponses) zeros(numel(binss{2}),8,4,nResponses)};
        %     latencyHists = {zeros(numel(binss{1})+1,8,4,nResponses) zeros(numel(binss{2})+1,8,4,nResponses)};
            isihs = cell(8,4,nResponses,2);

            for ii = 1:8
                for jj = 1:4
                    for kk = 1:nResponses
                        for ll = 1:2
                            isihs{ii,jj,kk,ll} = sparse([]);
                        end
                    end
                end
            end

            ns = zeros(8,4);

            
            randomIndices = randperm(nStimuli);
            trainIndices = sortedStimulusIndices(randomIndices(1:ceil(numel(randomIndices)/2)));
            testIndices = sortedStimulusIndices(randomIndices(ceil(numel(randomIndices)/2)+1:end));

            allCounts = zeros(numel(trainIndices),2,nResponses);
            allLatencies = zeros(numel(trainIndices),2,nResponses);
            allPhases = zeros(numel(trainIndices),1);
            allWidths = zeros(numel(trainIndices),1);

            nTrainTrials = numel(trainIndices);

            for ii = 1:nTrainTrials
                tic;
                stimulusIndex = gratings(trainIndices(ii));
                stimulusParams = getExtraParams{stimulusIndex};

                widthIndex = find(widths == stimulusParams.spatialPeriod*25/2);

                if isempty(widthIndex)
                    widthIndex = find(widths == stimulusParams.colour*25/2);
                end

                phaseIndex = find(phases == stimulusParams.phase*8);

                ns(phaseIndex,widthIndex) = ns(phaseIndex,widthIndex) + 1;
                allPhases(ii) = phaseIndex;
                allWidths(ii) = widthIndex;

                for jj = 0:1
                    interval = stimulusTimes(stimulusIndex+jj+[0 1]);

        %             tic;
        %             allSpikeIndices = find(sortedSpikeTimes >= interval(1) & sortedSpikeTimes < interval(2));
        %             toc;

                    for kk = 1:nResponses
        %                 cellSpikeIndices = allSpikeIndices( ...
        %                     sortedSpikeChannels(allSpikeIndices,1) == channel(1) & ...
        %                     sortedSpikeChannels(allSpikeIndices,2) == channel(2) & ...
        %                     sortedSpikeClusters(allSpikeIndices,1) == cluster);
        %                 tic;
                        cellSpikeIndices = find(cellSpikeTimes{kk} >= interval(1) & cellSpikeTimes{kk} < interval(2));


                        allCounts(ii,jj+1,kk) = numel(cellSpikeIndices);

                        if isempty(cellSpikeIndices)
                            allLatencies(ii,jj+1,kk) = Inf;
        %                     latencyHists{jj+1}(end,phaseIndex,widthIndex,kk) = ...
        %                         latencyHists{jj+1}(end,phaseIndex,widthIndex,kk) + 1;
                            continue;
                        end

                        spikeTimes = cellSpikeTimes{kk}(cellSpikeIndices);
                        spikeBins = floor((spikeTimes-interval(1))/binWidth)+1;

                        if strcmp(opts.latencytype,'burst')
                            if numel(spikeTimes) == 1
                                allLatencies(ii,jj+1,kk) = spikeTimes(1)-interval(1);
                            else
                                [~,~,burstStarts] = burstold(spikeTimes);

                                if isempty(burstStarts)
                                    allLatencies(ii,jj+1,kk) = spikeTimes(1)-interval(1);
                                else
                                    allLatencies(ii,jj+1,kk) = spikeTimes(burstStarts(1))-interval(1);
                                end
                            end
                        elseif strcmp(opts.latencytype,'peak')
                            psth = zeros(size(binss{jj+1}));
                            psth(spikeBins) = 1;
                            bw = 2*opts.latencyparam+1;
                            ifr = conv(psth,ones(bw,1)/bw,'same');
                            index = find(ifr == max(ifr),1);
                            allLatencies(ii,jj+1,kk) = index*binWidth;
                        else
                            allLatencies(ii,jj+1,kk) = spikeTimes(1)-interval(1);
                        end
        %                 latencyHists{jj+1}(index,phaseIndex,widthIndex,kk) = ...
        %                     latencyHists{jj+1}(index,phaseIndex,widthIndex,kk) + 1;

                        psths{jj+1}(spikeBins,phaseIndex,widthIndex,kk) = ...
                            psths{jj+1}(spikeBins,phaseIndex,widthIndex,kk) + 1;

        %                 tic;
        %                 previousSpike = find(sortedSpikeTimes < interval(1) & ...
        %                     sortedSpikeChannels(:,1) == channel(1) & ...
        %                     sortedSpikeChannels(:,2) == channel(2) & ...
        %                     sortedSpikeClusters == cluster,1,'last');
        %                 toc;

        %                 if isempty(previousSpike)
                        if cellSpikeIndices(1) == 1;
                            % assume an imaginary spike at the beginning of the
                            % recording to avoid having to deal with infinite ISIs
                            % or excluding trials from training: should be such a
                            % rare case as to not significantly decrease the
                            % performance of the algorithm
                            previousSpike = 0;
                        else
                            previousSpike = cellSpikeTimes{kk}(cellSpikeIndices(1)-1);
                        end

        %                 tic;
        %                 nextSpike = find(sortedSpikeTimes >= interval(2) & ...
        %                     sortedSpikeChannels(:,1) == channel(1) & ...
        %                     sortedSpikeChannels(:,2) == channel(2) & ...
        %                     sortedSpikeClusters == cluster,1);
        %                 toc;

        %                 if isempty(nextSpike)
                        if cellSpikeIndices(end) == numel(cellSpikeTimes{kk})
                            % assume an imaginary spike at the time of the last
                            % stimulus (end of recording would be more consistent
                            % but that isn't readily available without loading the
                            % MCD file again)
                            nextSpike = stimulusTimes(end);
                        else
                            nextSpike = cellSpikeTimes{kk}(cellSpikeIndices(end)+1);
                        end

                        intervals = diff([previousSpike; spikeTimes; nextSpike]);
                        intervalBins = floor(intervals/binWidth)+1;

                        isih = isihs{phaseIndex,widthIndex,kk,jj+1};

                        for ll = 1:numel(intervalBins)
                            intervalBin = intervalBins(ll);

                            if numel(isih) < intervalBin
                                isih(intervalBin) = 1;
                            else
                                isih(intervalBin) = isih(intervalBin) + 1;
                            end
                        end

                        isihs{phaseIndex,widthIndex,kk,jj+1} = isih;

        %                 toc;
                    end
                end
                fprintf('Evaluated trial %d of %d in %f seconds\n',ii,nTrainTrials,toc);
            end

            phasePDFs = zeros(8,4);

            maxCount = max(max(max(allCounts)));
            countPDFs = zeros(maxCount+1,2,nResponses);

            for ii = 1:nResponses
                for jj = 1:2
                    counts = allCounts(:,jj,ii);
                    countPDFs(:,jj,ii) = hist(counts,0:maxCount)/numel(counts);
                end
            end

            countPhaseConditionalPDFs = zeros(maxCount+1,8,nResponses,2,4);

            ratess = {zeros(size(psths{1})) zeros(size(psths{2}))};
            isirs = cell(size(isihs));

            kdePoints = 2^7;
            infLatencyPhaseConditionalPDFs = zeros(8,4,2,nResponses);
            latencyPhaseConditionalPDFSupports = zeros(kdePoints,8,4,2,nResponses);
            latencyPhaseConditionalPDFs = zeros(kdePoints,8,4,2,nResponses);

            for ii = 1:4
                for jj = 1:8
                    indices = allPhases == jj & allWidths == ii;
                    n = sum(indices);
                    phasePDFs(jj,ii) = n;

                    for kk = 1:2
                        bins = binss{kk};
                        psth = psths{kk};

                        t0 = 0;
                        t1 = bins(end)*binWidth;

                        for ll = 1:nResponses
                            tic;
                            counts = allCounts(indices,kk,ll);
                            countPhaseConditionalPDFs(:,jj,ll,kk,ii) = hist(counts,0:maxCount)/numel(counts);

                            latencies = allLatencies(indices,kk,ll);
                            finiteLatencies = latencies(isfinite(latencies));

                            infLatencyPhaseConditionalPDFs(jj,ii,kk,ll) = 1-numel(finiteLatencies)/numel(latencies);

                            if ~isempty(finiteLatencies)
                                [~,pdf,ts] = kde(latencies,kdePoints,t0,t1);
                                latencyPhaseConditionalPDFSupports(:,jj,ii,kk,ll) = ts;
                                latencyPhaseConditionalPDFs(:,jj,ii,kk,ll) = pdf/sum(pdf);
                            end

                            rate = (psth(:,jj,ii,ll)/n)/binWidth;
                            maxRate = max(rate);

                            if strcmpi(opts.ifrfit,'bars')
                                try
                                    fit = barsP(rate,[0 max(bins)],n,bp);
                                    rate = fit.mean;
                                catch err %#ok<NASGU>
                                    bw = opts.ifrparam;
                                    rate = conv(rate,normpdf(-(5*bw):(5*bw),0,bw),'same');
                                end
                            elseif strcmpi(opts.ifrfit,'kde')
                                tic;
                                PSTH = psth(:,jj,ii,ll);
                                data = zeros(sum(PSTH),1);
                                index = [0; cumsum(PSTH)];

                                for mm = 1:size(PSTH,1)
                                    data(index(mm)+1:index(mm+1)) = bins(mm)*ones(PSTH(mm),1);
                                end

                                [~,y,x] = kde(data,2^12,bins(1),bins(end));
                                rate = interp1q(x',y,bins);
                                rate(~isfinite(rate)) = 0;
                                rate = maxRate*rate/max(rate); % kde scales the data in an incomprehensible way, so this step is necessary to get it back within a sensible range
                                toc;
                            else
                                bw = opts.ifrparam;
                                rate = conv(rate,normpdf(-(5*bw):(5*bw),0,bw),'same');
                            end

                            ratess{kk}(:,jj,ii,ll) = rate;

                            isir = (isihs{jj,ii,ll,kk}/n)/binWidth;

                            if false %strcmpi(opts.ifrfit,'bars')
                                fit = barsP([(0:numel(isir)-1) isir],[0 numel(isir)-1],n);
                                isir = fit.mean(numel(isir)+1:end);
                            else
                                isir = conv(full(isir(:)),normpdf(-25:25,0,1),'same');
                            end

                            isirs{jj,ii,ll,kk} = isir;
                            toc;
                        end
                    end
                end
            end

            phasePDFs = phasePDFs./repmat(sum(phasePDFs),8,1);
%         else
%             load('lalala');
%         end
        
        logBinWidth = log(binWidth);
        
        logPhasePDFs = log(phasePDFs);
        
        logratess = cell(size(ratess));
        
        for ii = 1:numel(ratess)
            logratess{ii} = safelog(ratess{ii});
        end
        
        logisirs = cell(size(isirs));
        
        for ii = 1:numel(ratess)
            logisirs{ii} = safelog(isirs{ii});
        end

    %     save('posteriors','phasePDFs','-append');

    %     countPerformance = zeros(nResponses,4,2);
    %     latencyPerformance = zeros(nResponses,4,2);
    %     timingPerformance = zeros(nResponses,4,2);
    %     correlationPerformance = zeros(nResponses,4,2);
        nWidths = zeros(1,4);

        allTestWidths = zeros(numel(testIndices),1);
        allTestPhases = zeros(numel(testIndices),1);
        posteriors = zeros(8,numel(testIndices),nResponses,2,2);
        priors = zeros(8,numel(testIndices),nResponses,2,4);

        for ii = 1:numel(testIndices)
            stimulusIndex = gratings(testIndices(ii));
            stimulusParams = getExtraParams{stimulusIndex};

            widthIndex = find(widths == stimulusParams.spatialPeriod*25/2);

            if isempty(widthIndex)
                widthIndex = find(widths == stimulusParams.colour*25/2);
            end

            allTestWidths(ii) = widthIndex;

            phaseIndex = find(phases == stimulusParams.phase*8);

            allTestPhases(ii) = phaseIndex;

            nWidths(widthIndex) = nWidths(widthIndex) + 1;

            for jj = 0:1
                bins = binss{jj+1};
                rates = ratess{jj+1};

                interval = stimulusTimes(stimulusIndex+jj+[0 1]);
    %             allSpikeIndices = find(sortedSpikeTimes >= interval(1) & sortedSpikeTimes < interval(2));

                for kk = 1:nResponses
                    tic;
                    cellSpikeIndices = find(cellSpikeTimes{kk} >= interval(1) & cellSpikeTimes{kk} < interval(2));
    %                 cellSpikeIndices = allSpikeIndices( ...
    %                     sortedSpikeChannels(allSpikeIndices,1) == channel(1) & ...
    %                     sortedSpikeChannels(allSpikeIndices,2) == channel(2) & ...
    %                     sortedSpikeClusters(allSpikeIndices,1) == cluster);

                    fprintf('%d seconds spent finding spikes\n',toc);
                    tic;

                    % spike count code
                    count = numel(cellSpikeIndices);

                    pCount = (countPhaseConditionalPDFs(count+1,:,kk,jj+1,widthIndex))';
                    % disabled the denominator while debugging bedcause I
                    % forgot to save it - shouldn't make a difference
                    % anyway
                    posterior = safelog(pCount)+logPhasePDFs(:,widthIndex); %/countPDFs(count+1,jj+1,kk);

                    posteriors(:,ii,kk,jj+1,1) = posterior;

    %                 if find(posterior == max(posterior)) == phaseIndex
    %                     countPerformance(kk,widthIndex,jj+1) = countPerformance(kk,widthIndex,jj+1) + 1;
    %                 end

                    fprintf('%d seconds spent evaluating count code\n',toc);
                    tic;

                    % spike latency/timing/correlation code
                    pInf = infLatencyPhaseConditionalPDFs(:,widthIndex,jj+1,kk);
                    pLatency = zeros(8,1);

                    if isempty(cellSpikeIndices)
                        previousSpike = find(cellSpikeTimes{kk} < interval(1),1,'last');
                        previousSpike = cellSpikeTimes{kk}(previousSpike);
                        pLatency = pInf;
                    elseif cellSpikeIndices(1) == 1;
                        previousSpike = 0;
                    else
                        previousSpike = cellSpikeTimes{kk}(cellSpikeIndices(1)-1);
                    end

                    if isempty(previousSpike)
                        previousSpike = 0;
                    end

    %                 if cellSpikeIndices(end) == numel(cellSpikeTimes{kk})
    %                     nextSpike = stimulusTimes(end);
    %                 else
    %                     nextSpike = cellSpikeTimes{kk}(cellSpikeIndices(end)+1);
    %                 end

                    spikeTimes = [previousSpike; cellSpikeTimes{kk}(cellSpikeIndices)]-interval(1);

                    fprintf('%d seconds spent finding previous spike\n',toc);
                    tic;

                    if ~isempty(cellSpikeIndices)
                        spikeBins = floor(spikeTimes/binWidth)+1;
                        psth = zeros(size(bins));
                        psth(spikeBins(spikeBins > 0 & spikeBins <= numel(bins))) = 1;
                        ifr = conv(psth,ones(101,1)/101,'same');
                        index = find(ifr == max(ifr),1);
                        latency = index*binWidth;
                    end

                    isis = diff(spikeTimes);

                    fprintf('%d seconds spent finding response latencies and ISIs\n',toc);

                    timingPrior = zeros(8,1);
                    correlationPrior = zeros(8,1);

                    for ll = 1:8
    %                     tic;
                        if ~isempty(cellSpikeIndices)
                            ts = latencyPhaseConditionalPDFSupports(:,ll,widthIndex,jj+1,kk);
                            pdf = latencyPhaseConditionalPDFs(:,ll,widthIndex,jj+1,kk);
                            p = interp1q(ts,pdf,latency);

                            if isnan(p)
                                pLatency(ll) = 0;
                            else
                                pLatency(ll) = (1-pInf(ll))*p;
                            end
                        end

                        fprintf('%d seconds spent evaluating latency code\n',toc);
                        tic;

                        rate = rates(:,ll,widthIndex,kk);
                        exponent = -trapz(rate)*binWidth;
                        timingPrior(ll) = exp(exponent);

                        priors(ll,ii,kk,jj+1,1) = exponent;

                        % the transpose here is needed because
                        % spikeTimes(2:end) is an empty row matrix iff
                        % numel(spikeTimes) <= 1, but interp1q insists on a
                        % column matrix (even if that matrix is empty)
                        logifrs = interp1q(bins*binWidth, ...
                            logratess{jj+1}(:,ll,widthIndex,kk), ...
                            spikeTimes((2:end)'));
%                         timingPrior(ll) = timingPrior(ll)*prod(ifrs*binWidth);
                        
%                         tic;
                        priors(ll,ii,kk,jj+1,2) = sum(logifrs+logBinWidth);
%                         toc;

                        fprintf('%d seconds spent evaluating timing code\n',toc);
                        tic;

                        isir = isirs{ll,widthIndex,kk,jj+1};

                        % plus one because anything within the first binWidth
                        % of the trial goes in the first bin, so if the first
                        % pre-trial spike is within one binWidth of the trial
                        % start then it goes in the bin before the start of the
                        % trial
                        tau = bins-floor(spikeTimes(1)/binWidth)+1;

                        for mm = 2:numel(spikeTimes);
                            spikeTime = spikeTimes(mm);
                            spikeBin = floor(spikeTime/binWidth)+2;
                            tau(spikeBin:end) = tau(spikeBin:end)-tau(spikeBin);
                        end

                        tau = tau+1; % +1 because bin 0 is index 1

                        % assume the firing rate probability as a function of
                        % isi is zero outside of the range of all isis
                        % encountered in the training dataset
                        v = zeros(size(tau));
                        v(tau <= numel(isir)) = isir(tau(tau <= numel(isir)));
                        exponent = -trapz(rate.*v)*binWidth;
                        correlationPrior(ll) = exp(exponent);

                        priors(ll,ii,kk,jj+1,3) = exponent;
                        
                        logisir = logisirs{ll,widthIndex,kk,jj+1};

                        if isempty(isis)
                            % empty product is one
                            continue;
                        elseif isempty(logisir)
                            v = realmin*ones(size(isis));
                        else
                            v = interp1q((0:numel(logisir)-1)'*binWidth,logisir,isis);
                        end
                        

%                         correlationPrior(ll) = correlationPrior(ll)*prod(ifrs.*v);
                        
%                         tic;
                        priors(ll,ii,kk,jj+1,4) = sum(logifrs+v);
%                         toc;

                        fprintf('%d seconds spent evaluating correlation code\n',toc);
    %                     tic;
                    end

                    latencyPosterior = safelog(pLatency)+logPhasePDFs(:,widthIndex);
                    posteriors(:,ii,kk,jj+1,2) = latencyPosterior;

    %                 if find(latencyPosterior == max(latencyPosterior)) == phaseIndex
    %                     latencyPerformance(kk,widthIndex,jj+1) = latencyPerformance(kk,widthIndex,jj+1) + 1;
    %                 end

    %                 timingPosterior = timingPrior.*phasePDFs(:,widthIndex);
    %                 posteriors(:,ii,kk,widthIndex,jj+1,3) = timingPosterior;

    %                 if find(timingPosterior == max(timingPosterior)) == phaseIndex
    %                     timingPerformance(kk,widthIndex,jj+1) = timingPerformance(kk,widthIndex,jj+1) + 1;
    %                 end

    %                 correlationPosterior = correlationPrior.*phasePDFs(:,widthIndex);
    %                 posteriors(:,ii,kk,widthIndex,jj+1,4) = correlationPosterior;

    %                 if find(correlationPosterior == max(correlationPosterior)) == phaseIndex
    %                     correlationPerformance(kk,widthIndex,jj+1) = correlationPerformance(kk,widthIndex,jj+1) + 1;
    %                 end

    %                 fprintf('%d seconds spent predicting stimuli\n',toc);
                end
            end
        end
        
        save([fileDir '\' recording.dataFile '_bayes_' opts.ifrfit '_' num2str(opts.ifrparam) '_' opts.latencytype '_' num2str(opts.latencyparam) '.mat'],'priors','posteriors','allTestWidths','allTestPhases','nWidths','phasePDFs','logPhasePDFs','widths');
    else
%     
%     save('posteriors','priors','posteriors','allTestWidths','allTestPhases');
%     return;
    
        if ~iscell(recording)
            recordings = {recording};
            [allFileDir,allFilename] = getAnalysisOutputDir(recording);
        else
            recordings = recording;
            allFileDir = pwd;
            [~,allFilename] = fileparts(allFileDir);
        end
        
        allResponses = 0;
        for ii = 1:numel(recordings)
            recording = recordings{ii};
            
            [fileDir,filename,parentDir] = getAnalysisOutputDir(recording);
            
            posteriorsFile = [fileDir '\' filename '_bayes_' opts.ifrfit '_' num2str(opts.ifrparam) '_' opts.latencytype '_' num2str(opts.latencyparam) '.mat'];
            
            if ischar(parentDir)
                posteriorsFile = [parentDir '\' posteriorsFile]; %#ok<AGROW>
            end
            
            load(posteriorsFile);
            
            nResponses = size(posteriors,3); %#ok<NODEF>
            allResponses = allResponses + nResponses;
            
            if ii == 1
                allPriors = priors; %#ok<NODEF>
                allPosteriors = posteriors;
            else
                allPriors(:,:,end+(1:nResponses),:,:) = priors; %#ok<NODEF>
                allPosteriors(:,:,end+(1:nResponses),:,:) = posteriors;
            end
        end
        
        nResponses = allResponses;
        
        priors = allPriors;
        posteriors = allPosteriors;
        
        fileDir = allFileDir;
        filename = allFilename;
    end
    
    if all(logical(opts.posteriorsonly))
        return;
    end
    
%     nResponses = size(posteriors,3); %#ok<NODEF>
    
    if ndims(posteriors) == 6
        % when developing I accidentally saved a set of test data in the
        % wrong format; this corrects for that and should never happen in
        % practice
        posteriors = squeeze(max(posteriors,[],4));
        priors = squeeze(max(priors,[],4));
    end
    
%     nWidths = zeros(1,4);
    
    neuronDroppingStep = 5;
    
    if ceil(nResponses/neuronDroppingStep) == floor(nResponses/neuronDroppingStep)
        nNeurons = neuronDroppingStep:neuronDroppingStep:nResponses-neuronDroppingStep;
    else
        nNeurons = neuronDroppingStep:neuronDroppingStep:(neuronDroppingStep*floor(nResponses/neuronDroppingStep));
    end
    
    nRepeats = zeros(size(nNeurons));
    
    for ii = 1:numel(nNeurons)
        nRepeats(ii) = nchoosek(nResponses,nNeurons(ii));
    end
       
    maxRepeats = 10;
    nRepeats = min([nRepeats; maxRepeats*ones(size(nRepeats))]);
    
%     nNeurons = [1 nNeurons nResponses];
%     nRepeats = [nResponses nRepeats 1];
    
%     nNeurons = [1 floor(nResponses/2) nResponses];
%     nRepeats = [nResponses 1 1];

    nNeurons = nResponses;
    nRepeats = 1;
    
    performances = cell(size(nNeurons));
    
    for ii = 1:numel(nNeurons)
        performances{ii} = zeros(nRepeats(ii),4,2,4);
    end
    
    nTests = size(posteriors,2);
    
%     singleCellPerformance = zeros(nResponses,4,2,4);
%     allCellsPerformance = zeros(4,2,4);
    
    for ii = 1:nTests
        tic;
        widthIndex = allTestWidths(ii);
        phaseIndex = allTestPhases(ii);
        
%         nWidths(widthIndex) = nWidths(widthIndex) + 1;
        
        logPhasePDF = logPhasePDFs(:,widthIndex);
        
%         for jj = 1:nResponses
%             for kk = 1:2
%                 for ll = 1:2
%                     posterior = posteriors(:,ii,jj,kk,ll);
%                     
%                     if find(posterior == max(posterior)) == phaseIndex
%                         singleCellPerformance(jj,widthIndex,kk,ll) = ...
%                             singleCellPerformance(jj,widthIndex,kk,ll) + 1;
%                     end
%                     
%                     prior1 = priors(:,ii,jj,kk,2*ll-1);
%                     prior2 = priors(:,ii,jj,kk,2*ll);
%                     
%                     posterior = prior1.*prior2.*phasePDF;
%                     
%                     if find(posterior == max(posterior)) == phaseIndex
%                         singleCellPerformance(jj,widthIndex,kk,ll+2) = ...
%                             singleCellPerformance(jj,widthIndex,kk,ll+2) + 1;
%                     end
%                 end
%             end
%         end
    
        for mm = 1:numel(nNeurons)
            performance = performances{mm};
            n = nNeurons(mm);
            repeats = nRepeats(mm);
            
            for jj = 1:repeats
                if n == 1
                    indices = jj;
                elseif n == nResponses
                    indices = 1:nResponses;
                else
                    indices = randperm(nResponses);
                    indices = indices(1:n);
                end
                
                % TODO : this is a temporary fix, in reality it would be
                % better to recompute the posteriors using the new method
                allCellsPosterior = squeeze(sum(posteriors(:,ii,indices,:,:),3));
                allCellsPrior1 = squeeze(sum(priors(:,ii,indices,:,[1 3]),3));
                allCellsPrior2 = squeeze(sum(priors(:,ii,indices,:,[2 4]),3));
                
                for kk = 1:2
                    for ll = 1:2
                        posterior = allCellsPosterior(:,kk,ll);

                        % if multiple posteriors are the same (e.g. they're
                        % all zero), just pick the first one, as at least
                        % that way you might get it right by chance,
                        % whereas if you pick multiple ones the condition 
                        % will always evaluate false
                        if find(posterior == max(posterior),1) == phaseIndex
                            performance(jj,widthIndex,kk,ll) = ...
                                performance(jj,widthIndex,kk,ll) + 1;
                        end

                        prior1 = allCellsPrior1(:,kk,ll);
                        prior2 = allCellsPrior2(:,kk,ll);

                        posterior = prior1+prior2+logPhasePDF;

                        if find(posterior == max(posterior),1) == phaseIndex
                            performance(jj,widthIndex,kk,ll+2) = ...
                                performance(jj,widthIndex,kk,ll+2) + 1;
                        end
                    end
                end
            end
            
            performances{mm} = performance;
        end
        toc;
    end
    
    for ii = 1:numel(nNeurons)
        performance = performances{ii};
        performanceSize = size(performance);
        performanceSize(2) = 1;
        performance = performance./repmat(nWidths,performanceSize);
        performances{ii} = performance;
    end
    
%     singleCellPerformance = singleCellPerformance./repmat(nWidths,[nResponses 1 2 4]);
%     allCellsPerformance = allCellsPerformance./repmat(nWidths',[1 2 4]);
    
%     return;
%     
%     countPerformance = countPerformance./repmat(nWidths,[nResponses 1 2]);
%     timingPerformance = timingPerformance./repmat(nWidths,[nResponses 1 2]);
%     correlationPerformance = correlationPerformance./repmat(nWidths,[nResponses 1 2]);
    
%     save([fileDir '\' recording.dataFile '_bayes.mat'],'countPerformance','timingPerformance','correlationPerformance');
    performanceFile = [fileDir '\' filename '_bayes_' opts.ifrfit '_' num2str(opts.ifrparam) '_' opts.latencytype '_' num2str(opts.latencyparam) '.mat'];
    
    if exist(performanceFile,'file')
        save(performanceFile,'performances');
    else
        save(performanceFile,'performances','-append');
    end
    
%     return;
    
%     figure;
%     set(gcf,'Position',[0 0 1920 1080]);
    
%     for ii = 1:2
%         for jj = 1:4
%             subplot(2,4,4*(ii-1)+jj);
%             hist([countPerformance(:,jj,ii) timingPerformance(:,jj,ii) correlationPerformance(:,jj,ii)]);
%             hold on;
%             line([1/8 1/8],ylim,'Color','r');
%             
%             if ii == 2 && jj == 4
%                 legend({'Count' 'Timing' 'Correlation' 'Chance'});
%             end
%         end
%     end
%     
%     saveas(gcf,[fileDir '\' recording.dataFile '_coding strategies.png']);

%     load('performances');
    pairs = [1 2; 1 3; 3 4];
    codes = {'Count' 'Latency' 'Timing' 'Correlation'};
    figure;
    set(gcf,'Position',[0 0 1920 1080]);
    
    performance = performances{1};
    
    for hh = 1:size(pairs,1)
        clf;
        for ii = 1:2
            for jj = 1:4
                subplot(2,4,4*(ii-1)+jj);
                scatter(performance(:,jj,ii,pairs(hh,1)),performance(:,jj,ii,pairs(hh,2)));
%                 scatter(countPerformance(:,jj,ii),timingPerformance(:,jj,ii));
                hold on;
                fplot(@(x) x,xlim);
                line(xlim,[0.125 0.125],'LineStyle','--','Color','k');
                line([0.125 0.125],ylim,'LineStyle','--','Color','k');

                if ii == 2 && jj == 1
%                     xlabel('Fraction correct: count code');
%                     ylabel('Fraction correct: timing code');
                    xlabel(sprintf('Fraction correct: %s code',codes{pairs(hh,1)}));
                    ylabel(sprintf('Fraction correct: %s code',codes{pairs(hh,2)}));
                end
            end
        end
        
        saveas(gcf,[fileDir '\' filename '_' codes{pairs(hh,1)} '_vs_' codes{pairs(hh,2)} '_' opts.ifrfit '_' num2str(opts.ifrparam) '_' opts.latencytype '_' num2str(opts.latencyparam) '.png']);
    end
    
%     nNeurons = [1 5 27];
    meanPerformances = zeros(numel(nNeurons),4,2,4);
    
    for ii = 1:numel(nNeurons)
        performance = performances{ii};
        meanPerformances(ii,:,:,:) = mean(performance,1);
    end
    
    figure
    set(gcf,'Position',[0 0 1920 1080]);
    
    for ii = 1:4
        for jj = 1:2
            subplot(2,4,4*(jj-1)+ii);
            plot(nNeurons,squeeze(meanPerformances(:,ii,jj,:)));
            
            if ii == 4 && jj == 1
                legend(codes);
            end
            
            if ii == 1 && jj == 2
                xlabel('# Neurons');
                ylabel('Performance');
            end
        end
    end
    
    saveas(gcf,[fileDir '\' filename '_performance_vs_neurons_' opts.ifrfit '_' num2str(opts.ifrparam) '_' opts.latencytype '_' num2str(opts.latencyparam) '.png']);
    
    performance = performances{end};
    
    lineStyles = '-:';
    colours = 'rgby';
    
    figure;
    hold on;
    
    for ii = 1:2
        for jj = 1:4
            plot(widths,squeeze(performance(1,:,ii,jj)),'LineStyle',lineStyles(ii),'Color',colours(jj));
        end
    end
    
    legend(codes,'Location','SouthEast');
    xlabel('Bar Width');
    ylabel('Performance');
    set(gca,'XTick',widths);
    
    saveas(gcf,[fileDir '\' filename '_performance_vs_bar_width_new_' opts.ifrfit '_' num2str(opts.ifrparam) '_' opts.latencytype '_' num2str(opts.latencyparam) '.png']);
    
%     figure;
%     set(gcf,'Position',[0 0 1920 1080]);
%     
%     for ii = 1:2
%         for jj = 1:4
%             subplot(2,4,4*(ii-1)+jj);
%             scatter(countPerformance(:,jj,ii),correlationPerformance(:,jj,ii));
%             hold on;
%             fplot(@(x) x,xlim);
%             line(xlim,[0.125 0.125],'LineStyle','--','Color','k');
%             line([0.125 0.125],ylim,'LineStyle','--','Color','k');
%             
%             if ii == 2 && jj == 1
%                 xlabel('Fraction correct: count code');
%                 ylabel('Fraction correct: correlation code');
%             end
%         end
%     end
%     
%     saveas(gcf,[fileDir '\' recording.dataFile '_count_vs_correlation.png']);
%     
%     figure;
%     set(gcf,'Position',[0 0 1920 1080]);
%     
%     for ii = 1:2
%         for jj = 1:4
%             subplot(2,4,4*(ii-1)+jj);
%             scatter(timingPerformance(:,jj,ii),correlationPerformance(:,jj,ii));
%             hold on;
%             fplot(@(x) x,xlim);
%             line(xlim,[0.125 0.125],'LineStyle','--','Color','k');
%             line([0.125 0.125],ylim,'LineStyle','--','Color','k');
%             
%             if ii == 2 && jj == 1
%                 xlabel('Fraction correct: timing code');
%                 ylabel('Fraction correct: correlation code');
%             end
%         end
%     end
%     
%     saveas(gcf,[fileDir '\' recording.dataFile '_timing_vs_correlation.png']);
    return;
    
    for ii = 1:nResponses %#ok<UNRCH>
        cellRasters = squeeze(sortedRasters(ii,:,:));
        
        for jj = 1:4
            widthRasters = squeeze(cellRasters(jj,:));
            allRaster = cell(1,0);
            allPhases = [];
            allOnCounts = [];
            allOffCounts = [];
%             allStimulusLengths = zeros(8,2);
            
            for kk = 1:8
                widthRaster = widthRasters{kk}(:);
                allRaster = [allRaster; widthRaster];
                allPhases = [allPhases; kk*ones(size(widthRaster))];
                
                stimulusTiming = sortedTimings{jj,kk};
                
%                 for ll = 1:2
%                     allStimulusLengths(kk,ll) = max(diff(stimulusTiming(:,ll,:)));
%                 end
                
                for ll = 1:numel(widthRaster)
                    
                    row = widthRaster{ll};
                    
                    for mm = 1:2
                        count = numel(row(row > stimulusTiming(1,mm,ll) - stimulusTiming(1,1,ll) & row <= stimulusTiming(2,mm,ll) - stimulusTiming(1,1,ll)));
                        
                        if mm == 1
                            allOnCounts = [allOnCounts; count];
                        else
                            allOffCounts = [allOffCounts; count];
                        end
                    end
                end
            end
            
            indices = randperm(numel(allRaster));
            
            % spike count code
            countss = {allOnCounts allOffCounts};
            
            for kk = 1:2
                counts = countss{kk};
                trainCounts = counts(trainIndices);
                trainPhases = allPhases(trainIndices);
                pCount = hist(trainCounts,0:max(counts))/numel(trainCounts);
                pPhase = hist(trainPhases,1:8)/numel(trainPhases);
                
                cpdf = zeros(numel(0:max(counts)),8);
                
                for ll = 1:numel(trainIndices)
                    count = trainCounts(ll) + 1;
                    phase = trainPhases(ll);
                    
                    cpdf(count,phase) = cpdf(count,phase) + 1;
                end
                
                for ll = 1:8
                    cpdf(:,ll) = cpdf(:,ll)/sum(cpdf(:,ll));
                end
                
                performance = 0;
                
                for ll = 1:numel(testIndices)
                    index = testIndices(ll);
                    count = counts(index) + 1;
                    phase = allPhases(index);
                    
                    posterior = cpdf(count,:).*pPhase/pCount(count);
                    
                    if find(posterior == max(posterior)) == phase
                        performance = performance + 1;
                    end
                end
                
                countPerformance(ii,jj,kk) = performance/numel(testIndices);
            end

            % spike timing code
            xs = {(0:ceil(maxStimulusLengths(1)/binWidth))' (0:ceil(maxStimulusLengths(2)/binWidth))'};
            ys = {zeros(numel(xs{1}),8) zeros(numel(xs{2}),8)};
            ts = {(0:ceil((sum(maxStimulusLengths)+2)/binWidth))' (0:ceil((sum(maxStimulusLengths)+2)/binWidth))'};
            zs = {zeros(numel(ts{1}),8) zeros(numel(ts{2}),8)};
            ns = zeros(1,8);
            ms = zeros(2,8);
            
            trainPhases = allPhases(trainIndices);
            trainRaster = allRaster(trainIndices);
            correlationTrainingTrials = numel(trainIndices);
            
            for ll = 1:numel(trainIndices)
                row = trainRaster{ll};
                stimulusLengths = allStimulusLengths(trainIndices(ll),:);
                phase = trainPhases(ll);
                ns(phase) = ns(phase)+1;
                
                for mm = 1:2
                    spikeIndices = find(row > sum(stimulusLengths(1:mm-1)) & row <= sum(stimulusLengths(1:mm)));
                    spikes = row(spikeIndices) - sum(stimulusLengths(1:mm-1));
                    
                    spikeBins = floor(spikes/binWidth)+1;
                    
                    ys{mm}(spikeBins,phase) = ys{mm}(spikeBins,phase) + 1;
                    
                    precedingSpikes = spikeIndices - 1;
                    
                    if any(precedingSpikes < 1)
                        continue;
                    end
                    
                    intervals = spikes - (row(precedingSpikes) - sum(stimulusLengths(1:mm-1)));
                    
                    intervals = floor(intervals(isfinite(intervals))/binWidth)+1;
                    
                    zs{mm}(intervals,phase) = zs{mm}(intervals,phase) + 1;
                    
                    ms(mm,phase) = ms(mm,phase) + 1;
                end
            end
            
            vs = {zeros(size(ys{1})) zeros(size(ys{2}))};
            ws = {zeros(size(zs{1})) zeros(size(zs{2}))};
            
            for ll = 1:2
                x = xs{ll};
                t = ts{ll};
                y = ys{ll};
                y = (y./repmat(ns,size(y,1),1))/binWidth;
                z = zs{ll};
                z = (z./repmat(ms(ll,:),size(z,1),1))/binWidth;
                
                for mm = 1:8
%                     fit = barsP([x y(:,mm)],[0 max(x)],ns(mm));
%                     vs{ll}(:,mm) = fit.mean(numel(x)+1:end);
                    v = conv(y(:,mm),normpdf(-5:5,0,1),'same');
                    vs{ll}(:,mm) = v;
                    
%                     fit = barsP([t z(:,mm)],[0 max(t)],ns(mm));
%                     ws{ll}(:,mm) = fit.mean(numel(t)+1:end);
                    w = conv(z(:,mm),normpdf(-5:5,0,1),'same');
                    ws{ll}(:,mm) = w;
                end
            end    
            
            testPhases = allPhases(testIndices);
            testRaster = allRaster(testIndices);
            
            for ll = 1:2
                x = xs{ll};
                t = ts{ll};
                v = vs{ll};
                w = ws{ll};
                performance = 0;
                corrPerformance = 0;
                correlationCodeAttempts = numel(testIndices);
                
                for mm = 1:numel(testIndices)
                    phase = testPhases(mm);
                    pPhase = ns'/sum(ns);
                    
                    row = testRaster{mm};
                    stimulusLengths = allStimulusLengths(testIndices(ll),:);
                    
                    spikeIndices = find(row > sum(stimulusLengths(1:ll-1)) & row <= sum(stimulusLengths(1:ll)));
                    
                    if isempty(spikeIndices)
                        posterior = pPhase;
                        
                        for nn = 1:8
                            posterior = posterior*exp(-trapz(v(:,nn))*binWidth);
                        end
                        
                        if any(row <= sum(stimulusLengths(1:ll-1)))
                            lastSpike = row(find(row <= sum(stimulusLengths(1:ll-1)),1,'last')) - sum(stimulusLengths(1:ll-1));
                            
                            correlationPosterior = posterior*exp(-trapz(v(:,nn).*w((1:size(v,1))-floor(lastSpike/binWidth),nn))*binWidth);
                        else
                            % no idea if this is right, but in the case
                            % where we have no spikes and we have don't
                            % know when the last spike was, it seems
                            % intuitive that we'd have no more information
                            % than the timing code
                            correlationPosterior = posterior;
                        end
                        
                        doCorrelationCode = true;
                    else
                        if spikeIndices(1) < 2 || spikeIndices(end) >= numel(row)
                            doCorrelationCode = false;
                            spikes = [NaN; row(spikeIndices); NaN];
                        else
                            doCorrelationCode = true;
                            spikes = row([spikeIndices(1)-1; spikeIndices; spikeIndices(end)+1]) - sum(stimulusLengths(1:ll-1));
                        end

                        prior = zeros(8,1);
                        correlationPrior = zeros(8,1);

                        for nn = 1:8
                            prior(nn) = exp(-trapz(v(:,nn))*binWidth);
                            
                            if ~doCorrelationCode
                                continue;
                            end

                            spikeBins = ceil(spikes/binWidth);

                            lastSpikeTime = zeros(size(x));
                            
                            for oo = 2:numel(spikeBins)
                                idx = max(spikeBins(oo-1),1):min(spikeBins(oo)-1,numel(lastSpikeTime));
                                lastSpikeTime(idx) = spikeBins(oo-1);
                            end
                            
                            if spikeBins(end) == numel(lastSpikeTime)
                                lastSpikeTime(end) = numel(lastSpikeTime);
                            end
                            
                            tau = (1:numel(lastSpikeTime))'-lastSpikeTime+1;
    
                            correlationPrior(nn) = exp(-trapz(v(:,nn).*w(tau,nn))*binWidth);
                        end

                        for nn = 1:8
                            ifrs = interp1q(x*binWidth,v(:,nn),spikes(2:end-1));
                            prior(nn) = prior(nn)*prod(ifrs*binWidth);
                            
                            if ~doCorrelationCode
                                continue;
                            end

                            intervals = diff(spikes(1:end-1));
                            ifrs2 = interp1q(t*binWidth,w(:,nn),intervals);

                            correlationPrior(nn) = correlationPrior(nn)*prod(ifrs.*ifrs2);
                        end

                        posterior = prior.*pPhase;
                        correlationPosterior = correlationPrior.*pPhase;
                    end
                    
                    if find(posterior == max(posterior)) == phase
                        performance = performance + 1;
                    end
                    
                    if ~doCorrelationCode
                        correlationCodeAttempts = correlationCodeAttempts - 1;
                    elseif find(correlationPosterior == max(correlationPosterior)) == phase
                        corrPerformance = corrPerformance + 1;
                    end
                end
                
                timingPerformance(ii,jj,ll) = performance/numel(testIndices);
                correlationPerformance(ii,jj,ll) = corrPerformance/correlationCodeAttempts;
            end     
        end
    end
end