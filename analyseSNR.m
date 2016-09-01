 function analyseSNR(recording, varargin)
    opt = getopt('channel chlabel fetchremote=''no'' overwrite=''no'' lfp=''no'' maxtime=60 samplerate=25000 progressbar=''no'' window=2^15 freqs=0:5:300 maxlfpfreq=200',varargin);
       
    if isstruct(recording)
        if ~isfield(recording,'spont')
            return;
        elseif isempty(recording.spont) || ~recording.spont
            return;
        end
        
        channelInfo = getMetaInfo(recording.index,'spont',true);

        channels = parseChannelOptions(recording.index, 'channel',opt.channel,'chlabel',opt.chlabel);
        % TODO : test
        chlabels = vertcat(channelInfo(:).label);
        sampleRate = recording.sampleRate;
    else
        channels = 1:60;
        chlabels = num2str(channelIndexToMCSChannelNumber(channels'));
        sampleRate = opt.samplerate;
    end

    [dirname,filename] = getAnalysisOutputDir(recording);

    savefile = [dirname '\SNR.mat'];

    if ~strcmpi(opt.overwrite,'yes') && exist(savefile,'file')
        return;
    end

    nSpikess = cell(60,1);
    spikeRates = cell(60,1);
    noiseSpikes = zeros(60,1);
    noiseRate = zeros(60,1);
    clusterAmp = cell(60,1);
    clusterSNR = cell(60,1);
    noiseAmp = cell(60,1);
    noiseSNR = zeros(60,1);
    baseSTD = zeros(60,1);
    baseAmp = zeros(60,1);
    frequencies = cell(60,1);
    spectra = cell(60,1);
    spectrogramTimes = cell(60,1);
    spectrogramFreqs = cell(60,1);
    spectrograms = cell(60,1);
    psds = cell(60,1);
    max_clusters = -Inf;

    for cc = 1:length(channels)
        channel = channels(cc);
        tic;

        chlabel = chlabels(channel,:);

        spikesuffix = [filename '_channel_' chlabel];
        allspikefile = [dirname '\' spikesuffix '_spikes.mat'];
        allmcdspikefile = [dirname '\' spikesuffix '_MCD_trimmed_spikes.mat'];
        
        spikeprefix = [dirname '\times_' spikesuffix];
        clusterspikemanual = [spikeprefix '_concat.mat'];
        clusterspikeauto = [spikeprefix '_spikes.mat'];
        clustermcdspikefile = [spikeprefix '_MCD_trimmed_spikes.mat'];
        clusterspikefile = NaN;

        if exist(clusterspikemanual,'file')
            clusterspikefile = clusterspikemanual;
        elseif exist(clusterspikeauto,'file')
            clusterspikefile = clusterspikeauto;
        elseif exist(clustermcdspikefile,'file')
            clusterspikefile = clustermcdspikefile;
        end
        
        if ~any(isnan(clusterspikefile)) && exist(clusterspikefile,'file')
            load(clusterspikefile, 'spikes', 'cluster_class');
        elseif exist(allspikefile,'file')
            load(allspikefile, 'index', 'spikes');
        elseif exist(allmcdspikefile,'file')
            load(allmcdspikefile, 'index', 'spikes');
        else
            warning('Unable for find spike file for channel %d (%d/%d) in file %s',channel,cc,numel(channels),filename);
            continue;
        end
        
        hasSpikes = exist('spikes','var') && ~isempty(spikes); %#ok<NODEF>
        
        if ~exist('cluster_class','var')
            cluster_class = [ones(numel(index),1) index'];
        end
           
        if hasSpikes
            clusters = 1:max(cluster_class(:,1));

            nSpikes = zeros(length(clusters),1);

            snrs = [];
            amps = {};

            for cl = clusters;
                W = spikes(cluster_class(:,1) == cl,:);
                nSpikes(cl) = size(W,1);
                amps = [amps; {max(W,[],2) - min(W,[],2)}]; %#ok<AGROW>
%                 snrs = [snrs; sunersnr(W)];
            end

            nSpikess{channel} = nSpikes;

            max_clusters = max(length(snrs), max_clusters);

            clusterAmp{channel} = amps;
            clusterSNR{channel} = snrs;
            W = spikes(cluster_class(:,1) == 0,:);
            noiseSpikes(channel) = length(W);
            noiseAmp{channel} = max(W,[],2) - min(W,[],2);
            noiseSNR(channel) = sunersnr(W);

            if exist('channelInfo','var')
                maxtime = 1000*channelInfo(channel).segments/sampleRate;
            else
                maxtime = 300; % TODO : get time span from MCD file
            end

            spikeRates{channel} = nSpikes/maxtime;
            noiseRate(channel) = noiseSpikes(channel)/maxtime;
            
            clear spikes;
            clear cluster_class;
        end
        
        if ~strcmpi(opt.lfp,'yes')
            continue;
        end
        
        if exist('channelInfo','var')
            maxsamples = channelInfo(channel).segments;
        else
            maxsamples = Inf;
        end

        if hasSpikes
            % TODO : option to interactively choose segment of data
            % to analyse?

            times = cluster_class(:,2);
            dTimes = diff([0; times; maxtime]);

            maxDTimes = find(dTimes == max(dTimes),1);
            baseIndexCount = sampleRate*dTimes(maxDTimes)/1000;

            if maxDTimes == 1 || maxDTimes == length(dTimes)
                baseIndexCount = baseIndexCount - 64;
            else
                baseIndexCount = baseIndexCount - 128;
            end

            if baseIndexCount > opt.maxtime*sampleRate
                baseIndexCount = opt.maxtime*sampleRate;
            end

            if maxDTimes == 1
                baseStartIndex = 1;
            else
                baseStartIndex = sampleRate*times(maxDTimes-1)/1000 + 64;
            end

            if baseStartIndex + baseIndexCount > maxsamples
                baseIndexCount = maxSamples - baseStartIndex;
            end 
        else
            baseIndexCount = opt.maxtime*sampleRate;
            baseStartIndex = randi(maxsamples-baseIndexCount);
        end
        
        specIndexCount = 2^nextpow2(baseIndexCount);
        
        while maxsamples-specIndexCount < 1
            specIndexCount = specIndexCount/2;
        end
            
        specStartIndex = randi(maxsamples-specIndexCount);

        extractedDataFile = [dirname '\' filename '_channel_' chlabel '.mat'];
        
        if exist(extractedDataFile,'file')
            load(extractedDataFile,'channelData');
            alldata = channelData{1}; %#ok<USENS>
            
            if ~isempty(alldata)
                basedata = alldata(baseStartIndex:baseStartIndex+baseIndexCount-1);

                baseSTD(channel) = std(basedata)*1e6;
                baseAmp(channel) = max(abs(basedata-mean(basedata)))*1e6;

                clear basedata;
            else
                baseSTD(channel) = NaN;
                baseAmp(channel) = NaN;
            end

            specdata = alldata(specStartIndex:specStartIndex+specIndexCount-1);
        else
            [error,file] = ns_OpenFile([filename '.mcd']);
            
            if error
                warning('Could not find extracted data nor load raw MCD file for channel %s (%i/%i) in file %s', chlabel, cc, length(channels), filename); %#ok<WNTAG>
                continue;
            end
            
            [baseError,baseCount,basedata] = ns_GetAnalogData(file,channel,baseStartIndex,baseIndexCount); %#ok<NASGU>
            [specError,specCount,specdata] = ns_GetAnalogData(file,channel,specStartIndex,specIndexCount);
            
            if baseError || specError || baseCount < baseIndexCount || specCount < specIndexCount
                warning('Loading requested data segment failed for channel %s (%i/%i) in file %s', chlabel, cc, length(channels), filename); %#ok<WNTAG>
                continue;
            end
        end
        
        [b,a] = ellip(2,0.1,40,(2*opt.maxlfpfreq/sampleRate));
        lfpdata = filtfilt(b,a,specdata);
        
        freqSpectrum = fft(lfpdata);
        powSpectrum = abs(freqSpectrum).^2;
        frequencies{channel} = sampleRate/2*linspace(0,1,length(specdata)/2+1);
        spectra{channel} = powSpectrum(1:length(specdata)/2+1);

        [S,F,T,P] = spectrogram(specdata,hamming(opt.window),opt.window/2,opt.freqs,sampleRate);
        spectrogramTimes{channel} = T;
        spectrogramFreqs{channel} = F;
        spectrograms{channel} = S;
        psds{channel} = P;
        
        clear specdata;
        clear lfpdata;
        clear alldata;

        channelTime = toc;
        fprintf('Finished analysing SNR for channel %s (%i/%i) in file %s in %f seconds\n', chlabel, cc, length(channels), dirname, channelTime);
    end

    SNR = nan(60,max_clusters);
    meanSNR = zeros(60,1);

    for channel = channels
        snrs = clusterSNR{channel};

        if isempty(snrs)
            continue;
        end

        meanSNR(channel) = mean(snrs);
        SNR(channel,1:length(snrs)) = snrs;
    end

    save(savefile,'-v7.3','baseAmp','baseSTD','clusterAmp','noiseAmp','clusterSNR','meanSNR','noiseSNR','SNR','nSpikess','spikeRates','noiseSpikes','noiseRate','frequencies','spectra','spectrogramTimes','spectrogramFreqs','spectrograms','psds');
end