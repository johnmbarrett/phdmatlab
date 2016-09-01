% TODO : refactor
function condorSTA(spikeFile,timingFile,stimulusFile,nFrames,skipStimuli,verbose)
    if nargin < 5 || nargin > 6
        disp('Usage: condorSTA <spike file> <timing file> <stimulus file> <no. frames> <no. stimuli to ignore> <verbose?>');
        return;
    end
    
    fprintf([ ...
        'Running condorSTA with the following parameters:\n\n' ...
        'spikeFile: %s\n' ...
        'timingFile: %s\n' ...
        'stimulusFile: %s\n' ...
        'nFrames: %s\n' ...
        'skipStimuli: %s\n' ...
        'verbose: %s\n\n' ...
        ],spikeFile,timingFile,stimulusFile,nFrames,skipStimuli,verbose);

    nFrames = str2double(nFrames);
    skipStimuli = str2double(skipStimuli);
    verbose = nargin == 6 && strcmp(verbose,'true');
        
    if ~exist(timingFile,'file')
        error 'Timing File could not be found';
    end
    
    load(timingFile);
    
    if ~exist(spikeFile,'file')
        error 'Spike File could not be found';
    end
    
    load(spikeFile,'cluster_class');

    nTriggers = numel(stimulusTimes);
    stimulusTimes = stimulusTimes(skipStimuli+1:end)-recordingStartTime; %#ok<NODEF>
    nStimuli = numel(stimulusTimes);
    
    fprintf('Found %d stimulusFrames\n',nStimuli);

    sortedSpikeTimes = cluster_class(:,2)/100; %#ok<NODEF>
    
    cells = unique(cluster_class(:,1));
    nCells = numel(cells);
    
    fprintf('There are %d cells on this channel\n',nCells);
    
    nSpikes = zeros(nStimuli,nCells);
    
    for ii = 1:nStimuli-1
        if verbose
            tic;
        end

        interval = stimulusTimes(ii+[0 1]);

        spikes = find(sortedSpikeTimes > interval(1) & sortedSpikeTimes <= interval(2));

        if isempty(spikes)
            continue;
        end

        for jj = 1:nCells
            nSpikes(ii,jj) = numel(find(cluster_class(spikes,1) == cells(jj)));
        end
        
        if verbose
            toc;
        end
    end
    
    if ~exist(stimulusFile,'file')
        error 'Stimulus File could not be found';
    end
    
    reader = VideoReader(stimulusFile);
    
    stas = zeros(reader.Height,reader.Width,nFrames,nCells);
    
    pixels = double(read(reader,[1 nFrames]));
    pixels = squeeze(mean(pixels,3));
    
    for ii = nFrames:nStimuli-1
        tic;
        spikes = nSpikes(ii,:);
        spikes = reshape(spikes,1,1,1,nCells);
        spikes = repmat(spikes,[size(pixels,1) size(pixels,2) nFrames 1]);
        
        stas = stas + spikes.*repmat(pixels,[1 1 1 nCells]);
        
        pixels(:,:,1:end-1) = pixels(:,:,2:end);
        
        try
            pixels(:,:,end) = squeeze(mean(double(read(reader,ii+1)),3));
        catch err
            if isfield(err,'identifier') && strcmp(err.identifier,'MATLAB:read:invalidFrame')
                warning(['Failed to read stimulus %d of file %s.\n' ...
                    '# Frames\t%d\n' ...
                    '# Triggers\t%d\n' ...
                    '# Skipped\t%d\n' ...
                    '# Stimuli\t%d\n'], ...
                    ii+1, ...
                    stimulusFile, ...
                    reader.NumberOfFrames, ...
                    nTriggers, ...
                    skipStimuli, ...
                    nStimuli);
                error(err);
            end
        end
        
        toc;
    end
    
%     stas = zeros(reader.Height,reader.Width,nFrames,nCells);
%     
%     for ii = nFrames:nStimuli-1
%         if verbose
%             tic;
%         end
%         
%         pixels = double(read(reader,ii-nFrames+1));
%         
%         if size(pixels,3) > 1
%             pixels = mean(pixels,3);
%         end
%         
%         pixels = repmat(pixels,[1 1 nFrames nCells]);
%         
% %         if ii >= nFrames
%             spikes = nSpikes(ii:-1:ii-nFrames+1,:);
% %         else
% %             spikes = zeros(nFrames,nCells);
% %             spikes(end-ii+1:end,:) = nSpikes(ii:-1:1,:);
% %         end
%         
%         spikes = reshape(spikes,1,1,nFrames,nCells);
%         spikes = repmat(spikes,[size(pixels,1) size(pixels,2) 1 1]);
%         stas = stas + spikes.*pixels;
%         
%         if verbose
%             toc;
%         end
%     end
    
    NSpikes = sum(nSpikes);
    
    stas = stas./repmat(reshape(NSpikes,1,1,1,nCells),[size(stas,1) size(stas,2) nFrames 1]);
    stas(isnan(stas)) = 0;
    
%     stas2 = stas2./repmat(reshape(NSpikes,1,1,1,nCells),[size(stas2,1) size(stas2,2) nFrames 1]);
%     stas2(isnan(stas2)) = 0;
    
    timingMethodFileSuffixes = {'vsync_before' 'vsync_after' 'phototrigger'};
    averagingMethodFileSuffixes = {'all_spikes' 'first_spike'};
    
    tokens = regexp(spikeFile,'([\w\s]+\\)?(times_)?([\w\s]+)_channel_([0-9]+)','tokens');
    fileDir = tokens{1}{1};
    isClustered = ~isempty(tokens{1}{2});
    dataFile = tokens{1}{3};
    channel = tokens{1}{4};
    
    for jj = 1:nCells
        sta = squeeze(stas(:,:,:,jj));
        minSTA = min(min(min(sta)));
        maxSTA = max(max(max(sta)));
        sta = 255*(sta-minSTA)/(maxSTA-minSTA);
        sta = reshape(sta,size(sta,1),size(sta,2),1,size(sta,3));

        cluster = num2str(cells(jj));
        infix = '';

        if isClustered
            infix = '_clustered';
        end

        imwrite(sta,[fileDir 'sta_ ' dataFile '_channel_' channel '_cluster_' cluster infix '_' timingMethodFileSuffixes{3} '_' averagingMethodFileSuffixes{1} '.gif'],'LoopCount',1,'DelayTime',0.1);
        
%         sta = squeeze(stas2(:,:,:,jj));
%         minSTA = min(min(min(sta)));
%         maxSTA = max(max(max(sta)));
%         sta = 255*(sta-minSTA)/(maxSTA-minSTA);
%         sta = reshape(sta,size(sta,1),size(sta,2),1,size(sta,3));
% 
%         cluster = num2str(cells(jj));
%         infix = '';
% 
%         if isClustered
%             infix = '_clustered';
%         end
% 
%         imwrite(sta,[fileDir 'sta2_ ' dataFile '_channel_' channel '_cluster_' cluster infix '_' timingMethodFileSuffixes{3} '_' averagingMethodFileSuffixes{1} '.gif'],'LoopCount',1,'DelayTime',0.1);
    end
end