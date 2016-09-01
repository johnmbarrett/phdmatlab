function responsiveCells = analyseULEDSquareResponses(stimRecording,spontRecording,paramFilePrefix,varargin)
    functions = {@(varargin) NaN, NaN, @singleLEDs,@squares,@squaresNew,@(varargin) squarePSRH(true,varargin{:}),@(varargin) squarePSRH(false,varargin{:}),@frequency,@multiplePulses,@(varargin) movingBars(true,varargin{:}),@(varargin) movingBars(false,varargin{:}),@circles,@rings,@shapes};
    optionNames = {'nop' 'skip' 'single' 'squares' 'squaresnew' 'squarepsrh' 'squaresta' 'freq' 'multi' 'barpsrh' 'barsta' 'circles' 'rings' 'shapes'};
    
    optionSpec = '';
    
    for ii = 1:numel(optionNames)
        optionSpec = sprintf('%s %s=''no''',optionSpec,optionNames{ii});
    end
    
    options = getopt(sprintf('%s responselen=0.5 responsesig=0.05 threshmethod=''poiss'' figs=''yes'' analysespont=''no'' isistats=''no'' forceclustered=false blankttls=''yes'' allcells=''no'' responsivecells=NaN overwrite=''no'' functions=[] prefixes=[] suffixes=[]',optionSpec),varargin{:});
    
    if strcmpi(options.analysespont,'yes')
        [spontDir,spontFile] = getAnalysisOutputDir(spontRecording);
    
        try
            load(sprintf('%s\\SNR.mat',spontDir),'spikeRates');

            for ii = 1:numel(spikeRates) %#ok<NODEF>
                spikeRates{ii} = spikeRates{ii}*1000; %#ok<AGROW>
            end
        catch err %#ok<NASGU>
            % TODO : specify interval or read from file
            spikeRates = analyseSpont(spontFile,[0 300]);
        end
    else
        spikeRates = NaN;
    end
    
    [stimDir,stimFile] = getAnalysisOutputDir(stimRecording);
    
    load(sprintf('%s\\%s_vsync_times.mat',stimDir,stimFile),'vsyncTimes','recordingStartTime');
%     load(sprintf('%s\\%s_uled_timings.mat',stimDir,stimFile));
    
    ttls = vsyncTimes-recordingStartTime;
%     ttls = stimulusTimes-recordingStartTime;

    if strcmpi(options.isistats,'yes')
        firingStatsArgs = varargin;
        overwrite = find(strcmpi('overwrite',firingStatsArgs(1:2:end-1)));

        if ~isempty(overwrite)
            firingStatsArgs{overwrite*2} = 'no';
        end

        lambdas = getFiringStatistics(spontFile,firingStatsArgs{:});
    else
        lambdas = NaN;
    end
%     samples = [];
%     nCells = 0;
%     return;
    
    if ~isempty(options.functions) && numel(options.functions) == numel(options.suffixes)
        suffixes = options.suffixes;
        fns = cell(size(options.functions));
        
        for ii = 1:numel(options.functions)
            fns{ii} = functions{strcmp(options.functions{ii},optionNames)};
        end
    else
        fns = {};
        for ii = 1:numel(optionNames)
            option = options.(optionNames{ii});

            if isnumeric(option) && all(isfinite(option) & option >= 1 & option <= numel(functions)+1)
                for jj = 1:numel(option)
                    fns{option(jj)} = functions{ii}; %#ok<AGROW>
                end
            end
        end
        
        suffixes = 1:numel(fns);
    end
        
    if isempty(options.prefixes) || numel(options.prefixes) ~= numel(fns)
        prefixes = repmat({paramFilePrefix},size(options.functions));
    else
        prefixes = options.prefixes;
    end
    
    [~,currentDir] = fileparts(pwd);
    allResponsesFile = sprintf('%s_responsive_cells.mat',currentDir);
    newResponsesFile = sprintf('%s\\%s_uled_square_responses_newmethod.mat',stimDir,stimFile);
    someResponsesFile = sprintf('%s\\%s_uled_square_responses.mat',stimDir,stimFile);
    
    responsiveCells = [];
    if ~strcmp(options.allcells,'yes')
        if ischar(options.responsivecells) && exist(options.responsivecells,'file')
            load(options.responsivecells,'responsiveCells')
        elseif isnumeric(options.responsivecells) && ~any(isnan(reshape(options.responsivecells,numel(options.responsivecells),1)))
            responsiveCells = options.responsivecells;
        elseif exist(allResponsesFile,'file')
            load(allResponsesFile,'responsiveCells');
        elseif exist(newResponsesFile,'file')
            load(newResponsesFile,'responsiveCells');
        elseif exist(someResponsesFile,'file')
            load(someResponsesFile,'responsiveCells');
        end
            
        if ~isempty(responsiveCells) && ischar(responsiveCells)
            responsiveCells = [str2num(responsiveCells(:,1:2)) str2num(responsiveCells(:,4))]; %#ok<ST2NM>
        end
    end
    
    for ii = 1:numel(fns)
        if isnumeric(fns{ii}) && ~isempty(fns{ii}) && isnan(fns{ii})
            continue;
        end
        
        relativePath = sprintf('\\matlab\\stimuli\\%s_%d.txt',prefixes{ii},suffixes(ii));
        
        stimuli = fopen(sprintf('D:\\backup\\phd%s',relativePath),'r');
        
        if stimuli < 0
            stimuli = fopen(sprintf('V:\\retina\\John B\\phd backup%s',relativePath),'r');
            
            if stimuli < 0
                stimuli = fopen(sprintf('H:\\phd%s',relativePath),'r');
            end
        end
        
        [params,lastTTL] = parseStimFile(stimuli,options);

        if isa(fns{ii},'function_handle')
            returnVal = fns{ii}(stimFile,stimDir,params,ttls(1:lastTTL-1),spikeRates,lambdas,options,responsiveCells);
            
            if ~strcmp(options.allcells,'yes') && ~any(isnan(reshape(returnVal,numel(returnVal),1)))
                % TODO : this will stop working if any function other than
                % squares returns anything
                if  ~isempty(returnVal) && ischar(returnVal)
                    responsiveCells = [str2num(returnVal(:,1:2)) str2num(returnVal(:,4))]; %#ok<ST2NM>
                else
                    responsiveCells = returnVal;
                end
            end
        end
        
        if numel(ttls) >= lastTTL
            interval = ttls(lastTTL-[1 0]);

            if strcmp(options.analysespont,'yes')
                spikeRates = analyseSpont(spontFile,interval);
            end

            ttls = ttls(lastTTL:end);
        end
    end
end

function spikeRates = analyseSpont(recording,interval)
    t = diff(interval);
    spikeRates = cell(60,1);
    
    function fn(spikeTimes,channelIndex,~,cluster,varargin)
        spikes = find(spikeTimes > interval(1) & spikeTimes <= interval(2));
        spikeRates{channelIndex}(cluster) = numel(spikes)/t;
    end
    
    forEachChannel(recording,[],true,@fn);
end

function [X,Y,W,H,T,commandIndices,pws,xs,ys,widths] = parseSquaresFile(params)
    commandIndices = [params.r.index]';
    
    X = [params.r.X]';
    Y = [params.r.Y]';
    W = [params.r.W]';
    H = [params.r.H]';
    T = [params.r.T]';
    
    pws = unique(T);
    xs = unique(X);
    ys = unique(Y);
    widths = unique(W);
end

function [params,lastTTL] = parseStimFile(stimuli,options)
    params = struct( ...
        'b',struct('index',{},'T',{}), ...
        'e',struct('index',{},'X',{},'Y',{},'r',{},'R',{},'W',{},'T',{},'isContracting',{}), ...
        'g',struct('index',{},'F',{}), ...
        'i',struct('index',{},'T',{},'F',{}), ...
        'm',struct('index',{},'D',{},'T',{},'W',{}), ...
        'r',struct('index',{},'X',{},'Y',{},'W',{},'H',{},'T',{})  ...
        );
    
%     numericDataOffset = all(strcmp('g',commands) || strcmp('m',commands));
    
    lastTTL = 1;
    
    while ~feof(stimuli)
        s = fgetl(stimuli);
        
        if numel(s) < 1
            return;
        end
        
        command = s(1);
        
        isContracting = false;
        if strcmp(command,'c');
            isContracting = true;
            command = 'e';
        end
        
        params.(command)(end+1).index = lastTTL;
        
        switch command
            case 'b'
                if strcmp(options.blankttls,'yes')
                    A = sscanf(s,'b %d');
                    params.b(end).T = A(1);
                    lastTTL = lastTTL + 1;
                end
            case 'e'
                A = sscanf(s(2:end),' %d %d %d %d %d %d');
                params.e(end).X = A(1);
                params.e(end).Y = A(2);
                
                minR = A(3);
                params.e(end).r = minR;
                
                maxR = A(4);
                params.e(end).R = maxR;
                
                lastTTL = lastTTL + maxR - minR + 1;
                
                params.e(end).W = A(5);
                params.e(end).T = A(6);
                
                params.e(end).isContracting(end+1) = isContracting;
            case 'g'
                params.g(end).F = s(3:end);
                
                lastTTL = lastTTL + 1; % TODO : is this right?
            case 'i'
                tokens = regexp(s,'i\s+([0-9]+)\s+(.*)','tokens');
                params.i(end).T = str2double(tokens{1}{1});
                params.i(end).F = tokens{1}{2}; % no
                
                lastTTL = lastTTL + 1; % TODO : is this right?
            case 'm'
                A = sscanf(s,'m %s %d %d');
                params.m(end).D = char(A(1:end-2))';
                params.m(end).T = A(end-1);
                
                width = A(end);
                params.m(end).W = width;
                
                lastTTL = lastTTL + 16 - width;
            case 'r'
                A = sscanf(s,'r %d %d %d %d %d');
                params.r(end).X = A(1);
                params.r(end).Y = A(2);
                params.r(end).W = A(3);
                params.r(end).H = A(4);
                params.r(end).T = A(5);
                
                lastTTL = lastTTL + 1;
        end
    end
end
    
function nothing = singleLEDs(stimFile,stimDir,stimuli,ttls,spikeRates,varargin)
    nothing = NaN;
    [X,Y,~,~,T,commandIndices,pws] = parseSquaresFile(stimuli);
    
    figure;
    set(gcf,'Position',[0 0 1200 900],'PaperPositionMode','auto','Renderer','zbuffer');
            
    function fn(spikeTimes,channelIndex,channel,cluster,~)
        if numel(spikeRates{channelIndex}) < cluster || spikeRates{channelIndex}(cluster) == 0
            rate = numel(spikeTimes)/max(spikeTimes);
        else
            rate = spikeRates{channelIndex}(cluster);
        end

        rates = zeros(16,16,numel(pws));

        for jj = 1:numel(X)
            onset = ttls(commandIndices(jj));
            spikes = spikeTimes(spikeTimes > onset & spikeTimes <= onset+0.1);

            responseRate = numel(spikes)/0.1;
            rateChange = responseRate/rate - 1;

            x = X(jj)+1;
            y = Y(jj)+1;
            t = find(pws == T(jj));

            rates(y,x,t) = rates(y,x,t) + rateChange/10;
        end

        clf;
        [rows,cols] = subplots(numel(pws));

        maxRate = max(max(max(rates)));
        minRate = min(min(min(rates)));

        for jj = 1:numel(pws)
            subplot(rows,cols,jj);

            if minRate == 0
                C = rates(:,:,jj)/maxRate;
            else
                C = (rates(:,:,jj)-minRate)/(maxRate-minRate);
            end

            surf(squeeze(C));
            caxis([0 1]);
            view(2);
            shading interp;

            xlim([1 16]);
            ylim([1 16]);
            title(sprintf('PW = %d ms',pws(jj)));
        end

        filename = sprintf('%s\\%s_channel_%s_cluster_%d_single_uled_responses',stimDir,stimFile,channel,cluster);
        saveas(gcf,filename,'fig');
        saveas(gcf,filename,'png');
    end

    forEachChannel(stimFile,[],true,@fn);
end

function responsiveCells = squaresNew(stimFile,stimDir,stimuli,ttls,spikeRates,lambdas,options,responsiveCells)
    [X,Y,W,~,T,commandIndices,pws,xs,ys,widths] = parseSquaresFile(stimuli);
    nPWs = numel(pws);
    
    % TODO : this doesn't work in general
%     nConditions = numel(pws)*numel(xs)*numel(ys)*numel(widths);
    nStimuli = numel(commandIndices);
    groups = unique([X Y W],'rows');
    nConditions = size(unique([X Y W T],'rows'),1);
    nBlocks = nStimuli/nConditions;
    
    trainIndices = 1:2:nBlocks;
    testIndices = 2:2:nBlocks;
    nTests = numel(testIndices);
    
    maxPW = max(pws)/1000;
    
    stimulusTimes = ttls(commandIndices);
    preStimTimes = stimulusTimes-1;
    
    bws = [0.01 0.1];
    ts = [-0.1 0 0.1 0.2; -0.5 0 0.1 1.5] ;
    prefixes = {'raster' 'longraster'};
    
    function fn(spikeTimes,channelIndex,channel,cluster,~,spikes)
        if strcmp(options.figs,'yes')
            for ii = 1:2
                figs = rasterPlot(spikeTimes,stimulusTimes,T,ts(ii,:),[],true,bws(ii),[X Y W],{'PW' 'X' 'Y' 'Width'});

                for jj = 1:numel(figs)
                    fig = figs(jj);
                    set(fig,'Position',[0 0 1600 900]);
                    prefix = prefixes{ii};
                    filename = sprintf('%s\\%s_%s_channel_%s_cluster_%d_x_%d_y_%d_w_%d_uled_square_responses',stimDir,prefix,stimFile,channel,cluster,groups(jj,1),groups(jj,2),groups(jj,3));
                    saveas(fig,filename,'fig');
                    saveas(fig,filename,'png');
                    close(fig);
                end
            end
        end

        % TODO : don't need to recompute these every time, should really
        % save them all in a giant file somewhere to make reanalysis easier
        [~,~,~,spontSpikes] = rasterPlot(NaN,spikeTimes,preStimTimes,(1:nStimuli)',0,1,true,0.1,[],'PW',false);
        spontSpikes = spontSpikes(:);
        medianSpontSpikes = median(spontSpikes);
        meds = [meds; medianSpontSpikes];
        
        mu = mean(spontSpikes);
        mus = [mus; mu];
        
        % 0 std => Inf SNR, so put a lower bound on std equivalent to one
        % spontaneous spike in the entire recording
        sigma = max(std(spontSpikes),sqrt(1/numel(spontSpikes)));
        sigmas = [sigmas; sigma];
        
        [~,~,~,stimSpikes] = rasterPlot(NaN,spikeTimes,stimulusTimes,(1:nStimuli)',0,0.1,true,0.1,[],'PW',false);
        
        if reBootstrap
            sample = bootstrap(spontSpikes,numel(trainIndices),10000,true);
            maxStimSpikes = stimSpikes(T == pws(end));
            maxStimSpikes = maxStimSpikes(trainIndices);

            p = [p; sum(sample >= median(maxStimSpikes))/numel(sample)];
        end

        for ii = 1:numel(pws)
            pwSpikes = stimSpikes(T == pws(ii));
            pwSpikes = pwSpikes(testIndices);
            nSpikes(:,ii,end+(ii==1)) = pwSpikes;
            pSpikes(ii,end+(ii==1)) = sum(pwSpikes > medianSpontSpikes)/nTests;
            
            if sigma ~= 0
                % SNR = mean(S)/std(N), but don't have access to S, so
                % estimate as (mean(S+N)-mean(N))/std(N)
                snr(ii,end+(ii==1)) = (mean(pwSpikes)-mu)/sigma;
            else
                snr(ii,end+(ii==1)) = mean(pwSpikes);
            end
            
%             if mu ~= 0
%                 pSpikes(:,ii,end+(ii==1)) = (pwSpikes-mu)/mu;
%             else
%                 pSpikes(:,ii,end+(ii==1)) = pwSpikes;
%             end
%             
%             if sigma ~= 0
%                 zSpikes(:,ii,end+(ii==1)) = (pwSpikes-mu)/sigma;
%             else
%                 zSpikes(:,ii,end+(ii==1)) = pwSpikes;
%             end
%             
%             qSpikes(ii,end+(ii==1)) = sum(sample >= median(pwSpikes))/numel(sample);
        end

        channels = [channels; str2double(channel)];
        clusters = [clusters; cluster];
    end
    
    saveFile = sprintf('%s\\%s_uled_square_responses_newmethod.mat',stimDir,stimFile);
    
    reBootstrap = ~exist(saveFile,'file') || isempty(responsiveCells) || all(isnan(reshape(responsiveCells,numel(responsiveCells),1)));
    
    if strcmpi(options.overwrite,'yes') || ~exist(saveFile,'file')
        if reBootstrap
            p = [];
        end
        
        meds = [];
        mus = [];
        sigmas = [];
        nSpikes = zeros(nTests,nPWs,0);
        pSpikes = zeros(nPWs,0);
        snr = zeros(nPWs,0);
%         zSpikes = zeros(nTests,nPWs,0);
%         qSpikes = zeros(nPWs,0);
        channels = [];
        clusters = [];

        forEachChannel(stimFile,[],true,@fn,double(options.forceclustered));

        if reBootstrap
            responsive = fdrcorrect(p,0.05);

            responsiveCells = [channels(responsive) clusters(responsive)];
            nResponses = size(responsiveCells,1);
            responseIndices = 1:nResponses; 

            allNSpikes = nSpikes;
            allPSpikes = pSpikes;
            allSNR = snr;
            
            nSpikes = nSpikes(:,:,responsive);
            pSpikes = pSpikes(:,responsive);
            snr = snr(:,responsive);
%         zSpikes = zSpikes(:,:,responsive);
%         qSpikes = qSpikes(:,responsive);
        
        % largest p such that the probability of not observing an event
        % with probability p in n trials is 95%
%         qSpikes(qSpikes == 0) = 1-(1-0.05).^(1/numel(sample));

            save(saveFile,'channels','clusters','responsiveCells','allNSpikes','allPSpikes','allSNR','nSpikes','pSpikes','snr','meds','mus','sigmas'); %,'zSpikes','qSpikes');
        else
            load(saveFile,'channels','clusters');
            
            allNSpikes = nSpikes;
            allPSpikes = pSpikes;
            allSNR = snr;
            
            nResponses = size(responsiveCells,1);
            
            nSpikes = zeros(nTests,nPWs,nResponses);
            pSpikes = zeros(nPWs,nResponses);
            snr = zeros(nPWs,nResponses);
            
            cellIndices = [];
            responseIndices = [];
            
            for hh = 1:size(responsiveCells,1)
                index = find(ismember([channels clusters],responsiveCells(hh,:),'rows'));
                
                if ~isempty(index)
                    cellIndices = [cellIndices; index]; %#ok<AGROW>
                    responseIndices = [responseIndices; hh]; %#ok<AGROW>
                end
            end
            
            nSpikes(:,:,responseIndices) = allNSpikes(:,:,cellIndices);
            pSpikes(:,responseIndices) = allPSpikes(:,cellIndices);
            snr(:,responseIndices) = allSNR(:,cellIndices);
            
            save(saveFile,'channels','clusters','nSpikes','pSpikes','snr','allNSpikes','allPSpikes','allSNR','cellIndices','responseIndices','meds','mus','sigmas','-append'); %,'zSpikes','qSpikes');
        end
    else
        load(saveFile);
    end
    
    if strcmpi(options.figs,'yes')
        figure;
    end
    sigmoid = @(b,x) 1./(1+exp(-(x-b(1))/b(2)));
    beta0 = [0 1]; % pws(end)/2];
%     generalisedSigmoid = @(b,x) b(1)./(1+exp(-(x-b(2))/b(3)))+b(4);
%     gbeta0 = [1 0 1 0];
    
%     maxResponse = zeros(nResponses,1);
%     halfMaxPW = zeros(nResponses,1);
%     gain = zeros(nResponses,1);
%     backgroundFR = zeros(nResponses,1);


    thresholds = inf(nResponses,1);
    
    for hh = 1:nResponses
        responseIndex = responseIndices(hh);
        bp = nlinfit(pws,pSpikes(:,responseIndex),sigmoid,beta0);
        thresholds(responseIndex) = max(0,bp(1));
        
        if ~strcmpi(options.figs,'yes')
            continue;
        end
        
        clf;
%         hold on;
        
%         subplot(1,2,1);
%         plot(pws,qSpikes(:,hh));
%         
%         subplot(1,2,1);
%         hold on;
%         plot(pws,-log(qSpikes(:,hh)));
%         plot(pws,nSpikes(:,:,hh),'Color','k','LineStyle','none','Marker','o');
%         plot(pws,zSpikes(:,:,hh),'Color','r','LineStyle','none','Marker','o');
%         plot(pws,mean(pSpikes(:,:,hh)),'Color','c');
%         plot(pws,mean(zSpikes(:,:,hh)),'Color','m');
        
%         beta = nlinfit(pws,qSpikes(:,hh),sigmoid,beta0);
%         bn = nlinfit(kron(pws,ones(nTests,1)),reshape(nSpikes(:,:,hh),nTests*nPWs,1),generalisedSigmoid,gbeta0);
%         fplot(@(x) generalisedSigmoid(bn,x),xlim,'Color','b');
%         bp = nlinfit(kron(pws,ones(nTests,1)),reshape(pSpikes(:,:,hh),nTests*nPWs,1),sigmoid,beta0);
%         bz = nlinfit(kron(pws,ones(nTests,1)),reshape(zSpikes(:,:,hh),nTests*nPWs,1),sigmoid,beta0);
%         
%         subplot(1,2,2);
        hold on;
        
        plot(pws,pSpikes(:,hh),'Color','k','LineStyle','none','Marker','o');
        
        fplot(@(x) sigmoid(bp,x),xlim,'Color','b');
%         fplot(@(x) sigmoid(bz,x),xlim,'Color','m');
%         fplot(@(x) sigmoid(beta,x),xlim);
        
%         maxResponse(hh) = beta(1);
%         halfMaxPW(hh) = beta(2);
%         gain(hh) = 1./beta(3);
%         backgroundFR(hh) = beta(4);
%         thresholds(hh) = beta(2)-beta(3)*log(beta(A)/2-1);
        
        title(sprintf('Channel %d Cluster %d - Threshold: %5.2f ms',responsiveCells(hh,1),responsiveCells(hh,2),thresholds(hh)));
        xlabel('Pulse Width/ms');
        ylabel('Response Probability');
        ylim([0 1]);
        
        line(xlim,[0.5 0.5],'Color','k','LineStyle','--');
        line(bp([1 1]),[0 1],'Color','k','LineStyle','--');
        
        figFile = sprintf('%s\\thresholdsnew_%s_channel_%d_cluster_%d',stimDir,stimFile,responsiveCells(hh,1),responsiveCells(hh,2));
        saveas(gcf,figFile,'fig');
        saveas(gcf,figFile,'png');
    end
    
    save(saveFile,'thresholds','-append');
end

function responsiveCells = squares(stimFile,stimDir,stimuli,ttls,spikeRates,lambdas,options,responsiveCells)
    [X,Y,W,~,T,commandIndices,pws,xs,ys,widths] = parseSquaresFile(stimuli);
    
    % TODO : this doesn't work in general
%     nConditions = numel(pws)*numel(xs)*numel(ys)*numel(widths);
    groups = unique([X Y W],'rows');
    nConditions = size(unique([X Y W T],'rows'),1);
    nBlocks = numel(commandIndices)/nConditions;
    
    trainIndices = 1:2:nBlocks;
    testIndices = 2:2:nBlocks;
    
    maxPW = max(pws)/1000;
    
    plotFig = figure;
    figWidth = 100+numel(widths)*500;
    set(plotFig,'Position',[0 0 figWidth 550],'PaperPositionMode','auto','Renderer','zbuffer');
    
    threshFig = figure;
    
    responsiveCells = zeros(0,4);
    thresholdss = [];
    nSpikess = zeros(0,numel(pws));
    amplitudess = zeros(0,numel(pws));
    
    snrs = zeros(0,numel(pws));
%     snrs = zeros(nBlocks,numel(pws),0);
    allCells = zeros(0,2);
    
    bws = [0.01 0.1];
    ts = [-0.1 0 0.1 0.2; -0.5 0 0.5 1.5] ;
    prefixes = {'raster' 'longraster'};
    
    function fn(spikeTimes,channelIndex,channel,cluster,~,spikes)
        for ii = 1:2
            figs = rasterPlot(spikeTimes,ttls(commandIndices),T,ts(ii,:),[],true,bws(ii),[X Y W],{'PW' 'X' 'Y' 'Width'});
        
            for jj = 1:numel(figs)
                fig = figs(jj);
                set(fig,'Position',[0 0 1600 900]);
                prefix = prefixes{ii};
                filename = sprintf('%s\\%s_%s_channel_%s_cluster_%d_x_%d_y_%d_w_%d_uled_square_responses',stimDir,prefix,stimFile,channel,cluster,groups(jj,1),groups(jj,2),groups(jj,3));
                saveas(fig,filename,'fig');
                saveas(fig,filename,'png');
                close(fig);
            end
        end
        
        if any(isnan(spikeRates(:))) || numel(spikeRates{channelIndex}) < cluster || spikeRates{channelIndex}(cluster) == 0
            rate = numel(spikeTimes)/max(spikeTimes);
        else
            rate = spikeRates{channelIndex}(cluster);
        end
        if any(isnan(lambdas(:))) || numel(lambdas{channelIndex}) < cluster
            lambda = 0;
        else
            lambda = lambdas{channelIndex}(cluster);
        end
        
        ns = zeros(numel(ys),numel(xs),numel(pws),numel(widths));
        rates = zeros(nBlocks,numel(ys),numel(xs),numel(pws),numel(widths));
%         snr = zeros(nBlocks,numel(ys),numel(xs),numel(pws),numel(widths));
        responseProbs = zeros(numel(ys),numel(xs),numel(pws),numel(widths),2);
        thresholds = zeros(numel(ys),numel(xs),numel(widths));
        amplitudes = zeros(nBlocks/2,numel(ys),numel(xs),numel(pws),numel(widths));
        
        for ii = 1:numel(X)
            onset = ttls(commandIndices(ii));
            
%             spikesBefore = spikeTimes(spikeTimes > onset-options.responselen & spikeTimes <= onset);
            isInResponse = spikeTimes > onset & spikeTimes <= onset+options.responselen;
            spikesAfter = spikes(isInResponse,:);
            
            nSpikes = sum(isInResponse);

%             responseRate = numel(spikes)/maxPW;
%             rateChange = responseRate/rate - 1;

            x = find(xs == X(ii));
            y = find(ys == Y(ii));
            t = find(pws == T(ii));
            w = find(widths == W(ii));
            
            n = ns(y,x,t,w) + 1;
            ns(y,x,t,w) = n;
            
            isTest = ismember(n,testIndices);
            
            if isTest && ~isempty(spikesAfter)
                amplitudes(n/2,y,x,t,w) = median(max(spikesAfter,[],2)-min(spikesAfter,[],2));
            end
            
            rates(n,y,x,t,w) = nSpikes;
%             snr(n,y,x,t,w) = nSpikes/max(1,numel(spikesBefore));
            
            response = false;
            
            if nSpikes > 0
                if lambda == 0
                    response = true;
                else
                    response = (1-poisscdf(nSpikes-1,lambda*0.5)) < options.responsesig;
                end
            end
                       
            responseProbs(y,x,t,w,isTest+1) = responseProbs(y,x,t,w,isTest+1) + response;
        end
        
        snr = mean(rates)/sqrt(lambda*options.responselen);
%         snr = mean(rates)./std(rates);
%         snr(isnan(snr)) = 0;
        snrs(end+1,:) = snr(1,end,end,:,end); %#ok<SETNU>
%         snrs(:,:,end+1) = squeeze(snr(:,end,end,:,end)); %#ok<SETNU>
        allCells(end+1,:) = [str2double(channel) cluster]; %#ok<SETNU>
        
        responseProbs = responseProbs./(repmat(ns,[1 1 1 1 2])/2);
        
        medianRates = squeeze(median(rates(trainIndices,:,:,:,:),1));
        
        if (ndims(medianRates) == 5 && medianRates(1,1,end,end,1) < 1) ...
        || (ndims(medianRates) == 2 && medianRates(end,1) < 1) ... 
        || (all(squeeze(responseProbs(1,1,:,end,1)) < 0.5)) %#ok<ISMAT>
            % ignore all cells that fire on average less than one spike to
            % the strongest stimulus in the training data or which respond
            % less than half the time to all stimuli
            return;
        end
        
        medianRates = squeeze(median(rates(testIndices,:,:,:,:),1));
        
        responsiveCells = [responsiveCells; sprintf('%s %d',channel,cluster)];
        
        meanRates = squeeze(mean(rates(testIndices,:,:,:,:)));
%         ps = ones(size(meanRates)); %1-poisscdf(ceil(meanRates)-1,lambdas{channelIndex}(cluster));
        stdRates = 2*squeeze(std(rates(testIndices,:,:,:,:)))/sqrt(10);
        maxRate = max(max(max(max(meanRates+stdRates))));
        minRate = min(min(min(min(meanRates-stdRates))));
        minRate = min(minRate,0);
        rateRange = (maxRate-minRate)*1.05;
        
        if rateRange <= 0
            return;
        end
        
        figure(plotFig);
        clf;
        
        for ii = 1:numel(widths)
            validXs = unique(X(W == widths(ii)));
            validYs = unique(Y(W == widths(ii)));
            
            subplot('Position',[(55+525*(ii-1))/figWidth 0.10 475/figWidth 0.85]);
            hold on;
            
            xx = [0 numel(validXs)*(numel(pws)+1)];
            yy = [minRate minRate+rateRange*numel(validYs)];
            
            xlim(xx);
            ylim(yy);
            
            set(gca,'XTick',repmat(1:numel(pws),1,numel(validXs))+kron((0:numel(validXs)-1)*(numel(pws)+1),ones(1,numel(pws))));
            set(gca,'XTickLabel',repmat(pws(:),numel(validXs),1));
            
            tickStepChangePoints = repmat([1 2 5],1,7).*kron([0.1 1 10 100 100 1000 10000],ones(1,3));
            tickStep = tickStepChangePoints(find(tickStepChangePoints > rateRange,1))/10;
            
            ytick = ceil(minRate/tickStep)*tickStep:tickStep:floor(maxRate/tickStep)*tickStep;
            set(gca,'YTick',repmat(ytick,1,numel(validYs))+kron((0:numel(validYs)-1)*rateRange,ones(size(ytick))));
            set(gca,'YTickLabel',repmat(ytick,1,numel(validYs)));
            
            line(repmat((numel(pws)+1)*(1:numel(validXs)-1),2,1),repmat(yy',1,numel(validXs)-1),'Color','k');
            line(repmat(xx',1,numel(validYs)-1),repmat(minRate+rateRange*(1:numel(validYs)-1),2,1),'Color','k');
            
            if minRate < 0
                line(repmat(xx',1,numel(validYs)),repmat(rateRange*(0:numel(validYs)-1),2,1),'Color','k','LineStyle','--');
            end
            
            if rate > 0
                line(repmat(xx',1,numel(validYs)),repmat(rate*maxPW+rateRange*(0:numel(validYs)-1),2,1),'Color','k','LineStyle',':');
            end
            
            xlabel('Pulse width/ms');
            
            if ii == 1
                ylabel('# Spikes');
            end
            
            title(sprintf('Width: %dpx',widths(ii)));
            
            for jj = 1:numel(validXs)
                for kk = 1:numel(validYs)
                    x = xs == validXs(jj);
                    y = ys == validYs(kk);
                    % necessary for the significance plot because errorbar 
                    % plotting three 0x1 arrays plots the unity line, but 
                    % errorbar plotting three 1x0 arrays plots nothing, as
                    % expected
                    
%                     sample = samples{channelIndex}(cluster+1);
                    if numel(widths)*numel(validXs)*numel(validYs) == 1
                        r = meanRates;
                        s = stdRates;
%                         p = ps;
                    else
                        r = squeeze(meanRates(y,x,:,ii))';
                        s = squeeze(stdRates(y,x,:,ii))';

%                         p = squeeze(ps(y,x,:,ii));
                    end
                    
%                     p = p < 0.05/nCells;

%                     p = zeros(size(r));
                    
%                     for ll = 1:numel(r)
%                         p(ll) = sum(sample >= r(ll))/numel(sample) <= 0.05;
%                     end
                    
%                     if any(p)
%                        if ~significant
%                            responsiveChannels(end+1,:) = channel; %#ok<AGROW,SETNU>
%                            responsiveClusters = [responsiveClusters; cluster]; %#ok<AGROW>
%                            thresholds(:,end+1) = inf(numel(widths),1); %#ok<AGROW>
%                            significant = true;
%                        end
%                        
%                        thresholds(ii,end) = min(thresholds(ii,end),pws(find(p,1)));
%                     end
                    
                    figure(threshFig);
                    clf;
                    
                    if strcmpi(options.threshmethod,'poiss')
                        p = squeeze(responseProbs(y,x,:,ii,2));
                        
                        if isequal(p,zeros(size(p)))
                            threshold = Inf;
                        else
                            [expParams,expR] = nlinfit(pws,p,@(b,x) 1-exp(-(x-b(2))/b(1)),[1;0]);
                            [sigParams,sigR] = nlinfit(pws,p,@(b,x) 1./(1+exp(-(x-b(2))/b(1))),[1;0]);

                            if sum(expR.^2) < sum(sigR.^2)
                                a = expParams(2);
                                b = expParams(1);
                                threshold = a-b*log(0.5);
                            else
                                % logistic(0) = 0.5, so in this case the threshold is just the
                                % x offset
                                threshold = sigParams(2);
                            end
                        end
                        
                        if ii == numel(widths)
                            scatter(pws,squeeze(responseProbs(1,1,:,end,2)));

                            if isfinite(threshold)
                                hold on;
                                fplot(@(x) 1./(1+exp(-(x-sigParams(2))/sigParams(1))),[0 100],'Color','g');
                                fplot(@(x) 1-exp(-(x-expParams(2))/expParams(1)),[0 100],'Color','r');
                                line([threshold threshold],ylim,'Color','k','LineStyle','--');
                                line(xlim,[0.5 0.5],xlim,'Color','k','LineStyle','--');
                            end

                            xlabel('Pulse width/ms');
                            ylabel('Response probability');
                            ylim([0 1]);
                            title(sprintf('Channel %s cluster %d - threshold = %fms',channel,cluster,threshold));
                        end
                    elseif strcmpi(options.threshmethod,'hill')
                        if ismatrix(medianRates)
                            responses = medianRates;
                        else
                            responses = squeeze(medianRates(y,x,:,ii))';
                        end
                        
                        [threshold,exponent,maxR] = fitHillEquation(pws(:),responses(:));
                        
                        if ii == numel(widths)
                            hold on;
                            plot(pws,responses,'Color','b','LineStyle','none','Marker','o');
                        
                            if isfinite(threshold) && isfinite(exponent)
                                fplot(@(x) maxR*x.^exponent./(threshold^exponent + x.^exponent),[0 pws(end)],'Color','k');
                                line([threshold threshold],ylim,'Color','k','LineStyle','--');
                            end
                        end
                        
                        xlabel('Pulse width/ms');
                        ylabel('# Spikes');
                        title(sprintf('Channel %s cluster %d - threshold = %fms',channel,cluster,threshold));
                    elseif strcmpi(options.threshmethod,'sigmoid')
                        responses = squeeze(rates(testIndices,y,x,:,ii));
                        
                        sigmoid = @(b,x) b(1)./(1+exp(-(x-b(2))/b(3)))+b(4);
                        [params] = nlinfit( ...
                            repmat(pws,size(responses,1),1), ...
                            reshape(responses',numel(responses),1), ...
                            sigmoid, ...
                            [max(max(responses)) 0 1 min(min(responses))] ...
                            );
                        
                        threshold = params(2);
                        
                        if ii == numel(widths);
                            hold on;
                            plot(pws,responses,'Color','b','LineStyle','none','Marker','o');
                            fplot(@(x) sigmoid(params,x),[0 pws(end)],'Color','k');
                            line([threshold threshold],ylim,'Color','k','LineStyle','--');
                            xlabel('Pulse width/ms');
                            ylabel('# Spikes');
                            title(sprintf('Channel %s cluster %d - threshold = %fms',channel,cluster,threshold));
                        end
                    end
                    
                    thresholds(y,x,ii) = threshold;
                    
                    n = 1:numel(pws);
                    figure(plotFig);
                    errorbar(n+(jj-1)*(n(end)+1),r+rateRange*(kk-1),s);
                    
                    if any(pws > threshold)
                        errorbar(n(pws > threshold)+(jj-1)*(n(end)+1),r(pws > threshold)+rateRange*(kk-1),s(pws > threshold),'Color','r');
                    end

                    
                end
            end
        end
        
        threshold = thresholds(1,1,end);
        
        thresholdss(end+1) = max(0,threshold); %#ok<SETNU>
        nSpikess = [nSpikess; squeeze(rates(testIndices,1,1,:,end))];
        amplitudess = [amplitudess; squeeze(amplitudes(:,1,1,:,end))];
        
        filename = sprintf('%s\\thresholds_%s_channel_%s_cluster_%d_uled_square_responses',stimDir,stimFile,channel,cluster);
        saveas(threshFig,filename,'fig');
        saveas(threshFig,filename,'png');
        
        filename = sprintf('%s\\%s_channel_%s_cluster_%d_uled_square_responses',stimDir,stimFile,channel,cluster);
        saveas(plotFig,filename,'fig');
        saveas(plotFig,filename,'png');
    end

    forEachChannel(stimFile,[],true,@fn);
    
    save(sprintf('%s\\%s_uled_square_responses',stimDir,stimFile),'responsiveCells','thresholdss','nSpikess','amplitudess','pws','snrs','allCells');
end

function nothing = frequency(stimFile,stimDir,stimuli,ttls,spikeRates,lambdas,options,responsiveCells)
    nothing = NaN;
    [~,~,~,~,T,onsets,pws] = parseSquaresFile(stimuli);
    blanks = find([stimuli.b.T] == 5000);
    F = 1000./(T +  [stimuli.b(setdiff(1:numel(stimuli.b),blanks)).T]');
    freqs = unique(F);
    minPeriod = 1/max(freqs); 
    
    trainStarts = [1 blanks-(0:149)]';
    ts = T(trainStarts(1:end-1));
    fs = F(trainStarts(1:end-1));
    nTrains = numel(trainStarts)-1;
    stimulusTimes = ttls(onsets(trainStarts(1:end-1)));
    
    allCovNSpikes = [];
    allSpikesPerPulse = [];
    allFirstVLast = [];
    
    cells = zeros(0,2);
    
    function fn(spikeTimes,~,channel,cluster,~,~)
        if ~isempty(responsiveCells) && ~ismember([str2double(channel) cluster],responsiveCells,'rows')
            return;
        end
        
        if strcmp(options.figs,'yes') 
            [rasterFigs,~,~,~,~,subs] = rasterPlot(spikeTimes,stimulusTimes,fs,[-0.5 0 5 5.5],[],true,minPeriod,ts,{'Frequency','Pulse Width'});

            for ii = 1:numel(rasterFigs)
                set(rasterFigs(ii),'Position',[0 0 1600 900]);
            end
        end
        
        trials = zeros(numel(pws),numel(freqs));
        allNSpikes = zeros(size(onsets,1),4);
        allLatency = zeros(size(onsets,1),4);
        
        for ii = 1:nTrains
            n = trainStarts(ii);
            pw = T(n);
            freq = F(n);
            period = 1/freq;
            responseWindow = min(period,0.15);
            
            t = find(pw == pws);
            f = find(freq == freqs);

            trials(t,f) = trials(t,f) + 1;
            stimulusInfo = [t,f,trials(t,f)];
            
            pulses = onsets(n:trainStarts(ii+1)-1);
            
            assert(numel(pulses) == 5*freq);
            
            pulseTimes = ttls(pulses);
            
            for jj = 1:numel(pulseTimes)
                spikes = spikeTimes(spikeTimes > pulseTimes(jj) & spikeTimes <= pulseTimes(jj) + responseWindow);

                if isempty(spikes)
                    allNSpikes(n+jj-1,:) = [0 stimulusInfo];
                    allLatency(n+jj-1,:) = [Inf stimulusInfo];
                else
                    allNSpikes(n+jj-1,:) = [numel(spikes) stimulusInfo];
                    allLatency(n+jj-1,:) = [spikes(1) - pulseTimes(jj) stimulusInfo];
                end
            end
        end
        
        maxNTrials = max(max(trials));
        
        if isempty(allCovNSpikes)
            allCovNSpikes = nan(maxNTrials,numel(pws),numel(freqs),0);
            allSpikesPerPulse = nan(maxNTrials,numel(pws),numel(freqs),0);
            allFirstVLast = nan(maxNTrials,numel(pws),numel(freqs),0);
        end
        
        covNSpikes = nan(maxNTrials,numel(pws),numel(freqs));
        spikesPerPulse = nan(maxNTrials,numel(pws),numel(freqs));
        
        firstVLast = nan(maxNTrials,numel(pws),numel(freqs));
        
%         tau = nan(numel(pws),numel(freqs),3);
%         slope = nan(numel(pws),numel(freqs),3);
%         intercept = nan(numel(pws),numel(freqs),3);
        
%         spearRho = nan(numel(pws),numel(freqs));
%         spearPVal = nan(numel(pws),numel(freqs));
%         covLatency = zeros(maxNTrials,numel(pws),numel(freqs));
        
        for ii = 1:numel(freqs)
            freq = freqs(ii);
            
            pulsesPerTrain = 5*freq;
            
            for jj = 1:numel(pws)
                nTrials = trials(jj,ii);
                
                if nTrials == 0
                    continue;
                end
                
                nSpikes = zeros(pulsesPerTrain,nTrials);
                
                for kk = 1:nTrials;
                    index = allNSpikes(:,2) == jj & allNSpikes(:,3) == ii & allNSpikes(:,4) == kk;
                    
                    if sum(index) == 0
                        continue;
                    end
                    
                    nSpikes(:,kk) = allNSpikes(index,1);
                end
                
                mu = mean(nSpikes);
                
                spikesPerPulse(:,jj,ii) = mu;
                
                sigma = std(nSpikes);
                noSpikes = mu == 0;
                
                covNSpikes(noSpikes,jj,ii) = 0;
                covNSpikes(~noSpikes,jj,ii) = sigma(~noSpikes)./mu(~noSpikes);
                
                firstVLast(:,jj,ii) = diff(nSpikes([1 end],:));
                
%                 tdata = repmat((0:pulsesPerTrain-1)'*(1/freq),nTrials,1);
%                 xdata = reshape(nSpikes,pulsesPerTrain*nTrials,1);
%                 
%                 modelfun = @(b,x) b(2)*exp(b(1)*x)+b(3);
%                 [params,residuals,jacobian] = nlinfit(tdata,xdata,modelfun,[0; 0; 0]);
                
%                 if strcmp(options.figs,'yes') 
%                     axes(subs(ii,jj)); %#ok<LAXES>
%                     hold on;
%                     fplot(@(x) modelfun(params,x),xlim,'Color','r');
%                 end
%                 
%                 parci = nlparci(params,residuals,'jacobian',jacobian);
%                 
%                 growing = sign(params(1,1)) == sign(params(2,1));
%                 if growing == (params(1,1) >= 0)
%                     params = [params params-parci(:,1) parci(:,2)-params]; %#ok<AGROW>
%                 else
%                     params = [-params params-parci(:,2) parci(:,1)-params];
%                 end
%                 
%                 tau(jj,ii,:) = params(1,:);
%                 slope(jj,ii,:) = params(2,:);
%                 intercept(jj,ii,:) = params(3,:);
                
%                 if sigma == zeros(size(sigma))
%                     spearRho(jj,ii) = 0;
%                     spearPVal(jj,ii) = 1;
%                 else
%                     [rho,pval] = corr(repmat((1:pulsesPerTrain)',nTrials,1),reshape(nSpikes,pulsesPerTrain*nTrials,1),'type','Spearman');
%                     spearRho(jj,ii) = rho;
%                     spearPVal(jj,ii) = pval;
%                 end
                
                continue;
                
                ts = 0:(1/freq):5-(1/freq);
                
                mu = squeeze(mean(spikesPerPulse(:,:,jj)));
%                 sigma = squeeze(std(spikesPerPulse(:,:,jj))/sqrt(10));

                subplot(2,3,jj);
                hold on;
                plot(ts,mu,'Color',colours(ii,:));
                xlim([-1 5]);
                
                for kk = 1:pulsesPerTrain
                    latencies = pulseLatencies(:,kk,jj);
                    latencies = latencies(isfinite(latencies));
                    
                    if isempty(latencies)
                        mu(kk) = NaN;
%                         sigma(kk) = NaN;
                    else
                        mu(kk) = mean(latencies);
%                         sigma(kk) = std(latencies)/sqrt(numel(latencies));
                    end
                end

                subplot(2,3,3+jj);
                hold on;
                plot(ts,mu,'Color',colours(ii,:));
                xlim([-1 5]);
            end
            
%             meanNSpikes(:,ii) = mean(reshape(spikesPerPulse,10*pulsesPerTrain,numel(pws)));
%             stdNSpikes(:,ii) = std(reshape(spikesPerPulse,10*pulsesPerTrain,numel(pws)));
        end
        
        allCovNSpikes(:,:,:,end+1) = covNSpikes;
        allSpikesPerPulse(:,:,:,end+1) = spikesPerPulse;
        allFirstVLast(:,:,:,end+1) = firstVLast;
        
        cells(end+1,:) = [str2double(channel) cluster];
        
        if strcmp(options.figs,'yes')
            colours = distinguishable_colors(numel(pws)+3);
            figPrefix = sprintf('%s\\freqraster_%s_channel_%s_cluster_%d',stimDir,stimFile,channel,cluster);
            
            for ii = 1:numel(rasterFigs);
                fig = rasterFigs(ii);
                figFile = sprintf('%s_pw_%d',figPrefix,pws(ii));
                saveas(fig,figFile,'fig');
                saveas(fig,figFile,'png');
                close(fig);
            end
                    
            datas = {covNSpikes spikesPerPulse};
            prefixes = {'covvsfreq' 'spikesvsfreq'};
            ylabels = {'Coefficient of Variation' 'Spikes Per Pules'};
            
            for hh = 1:2
                fig = figure('Visible','off');
                hold on;

                data = datas{hh};
                for ii = 1:numel(pws)
                    mu = mean(squeeze(data(:,ii,:)));
                    sigma = std(squeeze(data(:,ii,:)));

                    errorbar(freqs,mu,sigma,'Color',colours(ii,:));
                end

                set(gca,'XScale','log','XTick',freqs);
                xlim([0.75 75]);
                xlabel('Frequency/Hz');
                ylabel(ylabels{hh});

                legend(cellstr([num2str(pws) repmat('ms',numel(pws),1)]),'Location','SouthEast');

                title(sprintf('Channel %s Cluster %d',channel,cluster));

                figFile = sprintf('%s\\%s_%s_channel_%s_cluster_%d',stimDir,prefixes{hh},stimFile,channel,cluster);
                saveas(fig,figFile,'fig');
                saveas(fig,figFile,'png');
                close(fig);
            end
            
            fig = figure; %('Visible','off');
            hold on;
            
            h = line([0 max(freqs)],[0 0],'Color','k','LineStyle','--');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
            
%             sig = spearPVal < 0.05;
%             verySig = spearPVal < 0.01;
%             veryVerySig = spearPVal < 0.001;
%             
%             hs = zeros(3,1);

%             tauM = tau(:,:,1);
%             
%             for ii = 1:numel(pws)
%                 plot(freqs,tauM(ii,:),'Color',colours(ii,:))
%             end
%             
%             yy = ylim;
%             clf;
%             hold on;
%             
%             tauL = tau(:,:,2);
%             tauU = tau(:,:,3);

            sig = false(numel(pws),numel(freqs));
            verySig = false(numel(pws),numel(freqs));
            veryVerySig = false(numel(pws),numel(freqs));
            
            for ii = 1:numel(pws)
                fvl = prctile(squeeze(firstVLast(:,ii,:)),[50 25 75]);
                fvl = [fvl; fvl(1,:)-fvl(2,:); fvl(3,:)-fvl(1,:)]; %#ok<AGROW>
                errorbar(freqs,fvl(1,:),fvl(2,:),fvl(3,:),'Color',colours(ii,:))
                
                for jj = 1:numel(freqs)
                    if ~any(isfinite(firstVLast(:,ii,jj)))
                        continue;
                    end
                    
                    p = signrank(firstVLast(:,ii,jj));
                    sig(ii,jj) = p < 0.05;
                    verySig(ii,jj) = p < 0.01;
                    veryVerySig(ii,jj) = p < 0.001;
                end
                
%                 errorbar(freqs,tauM(ii,:),tauL(ii,:),tauU(ii,:),'Color',colours(ii,:))
%                 plot(freqs,spearRho(ii,:),'Color',colours(ii,:));
%                 
                hs(1) = scatter(freqs(sig(ii,:)),fvl(1,sig(ii,:)'),'MarkerEdgeColor',colours(end-2,:));
                hs(2) = scatter(freqs(verySig(ii,:)),fvl(1,verySig(ii,:)'),'MarkerEdgeColor',colours(end-1,:));
                hs(3) = scatter(freqs(veryVerySig(ii,:)),fvl(1,veryVerySig(ii,:)'),'MarkerEdgeColor',colours(end,:));
                
                for jj = 1:3
                    set(get(get(hs(jj),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
            end
            
            xlabel('Frequency/Hz');
            ylabel('First Pulse Spikes - Last Pulse Spikes');
%             ylabel('Growth/Decay Rate');
%             ylabel('Spearman''s {\rho}');
            
            for ii = 2:-1:0
                scatter(-1,-1,'MarkerEdgeColor',colours(end-ii,:));
            end
            
            set(gca,'XScale','log','XTick',freqs);
            xlim([0.75 75]);
%             ylim(yy);
            
%             validUpper = find(tauU < 1e3 & isfinite(tauU));
%             validLower = find(tauL < 1e3 & isfinite(tauL));
%             [~,maxU] = max(tauU(validUpper));
%             [~,maxL] = max(tauL(validLower));
%             
%             ylim(1.1*[ ...
%                 tauM(validLower(maxU))-tauL(validLower(maxL)) ...
%                 tauM(validUpper(maxU))+tauU(validUpper(maxU)) ...
%                 ]);

%             legend(cellstr([num2str(pws) repmat('ms',numel(pws),1)])');
            legend([cellstr([num2str(pws) repmat('ms',numel(pws),1)])' {'p < 0.05' 'p < 0.01' 'p < 0.001'}],'Location','NorthEast');

            title(sprintf('Channel %s Cluster %d',channel,cluster));

            figFile = sprintf('%s\\firstvlast_%s_channel_%s_cluster_%d',stimDir,stimFile,channel,cluster);
            saveas(fig,figFile,'fig');
            saveas(fig,figFile,'png');
            close(fig);
        end
    end

    forEachChannel(stimFile,[],true,@fn);
    
    save(sprintf('%s\\%s_frequency_responses',stimFile,stimDir),'allCovNSpikes','allSpikesPerPulse','allFirstVLast','cells');
end

function nothing = multiplePulses(stimFile,stimDir,stimuli,ttls,spikeRates,lambdas,options,responsiveCells)
    nothing = NaN;
    
    stimulusEnds = find([stimuli.b.T] > 1000);
    stimulusStarts = [1 stimulusEnds(1:end-1)+1];
    stimulusTimes = ttls(stimulusStarts);
    nStimuli = numel(stimulusEnds);
    
    totalPW = zeros(nStimuli,1);
    innerPW = zeros(nStimuli,1);
    pausePW = zeros(nStimuli,1);
    
    for hh = 1:nStimuli
        pulses = stimuli.r(stimulusStarts(hh):stimulusEnds(hh));
        
        pulseTimes = [pulses.T];
        
        totalPW(hh) = sum(pulseTimes);
        
        assert(numel(unique(pulseTimes)) == 1);
        
        innerPW(hh) = pulseTimes(1);
        
        pauses = stimuli.b(stimulusStarts(hh):stimulusEnds(hh)-1);
        
        if isempty(pauses)
            pausePW(hh) = 0;
        else
            pauseTimes = [pauses.T];

            assert(numel(unique(pauseTimes)) == 1);

            pausePW(hh) = pauseTimes(1);
        end
    end
    
    relPW = innerPW./totalPW;
    relPWs = unique(relPW);
    pausePWs = unique(pausePW);
    totalPWs = unique(totalPW);
    
    [rows,cols] = subplots(numel(pausePWs));
    colours = distinguishable_colors(numel(relPWs)-1);
    
    allConditions = [totalPW relPW pausePW];
    conditions = unique(allConditions,'rows');
    nConditions = size(conditions,1);
    nTrials = ceil(nStimuli/nConditions);
    
    conditionIndices = zeros(nStimuli,1);
    
    for hh = 1:nStimuli
        conditionIndices(hh) = find(ismember(conditions,allConditions(hh,:),'rows'));
    end
    
    textConditions = cellstr([              ...
        repmat('T ',nConditions,1)          ...
        num2str(conditions(:,1),'%02d')     ...
        repmat('/I ',nConditions,1)         ...
        num2str(conditions(:,2),'%4.3f')    ...
        repmat('/P ',nConditions,1)         ...
        num2str(conditions(:,3),'%02d')     ...
        ]);
    
    allNSpikes = zeros(nTrials,nConditions,0);
    cells = zeros(0,2);
    
    function fn(spikeTimes,~,channel,cluster,~,~)
        if ~isempty(responsiveCells) && ~ismember([str2double(channel) cluster],responsiveCells,'rows')
            return;
        end
        
        if strcmp(options.figs,'yes') 
            rasterFig = rasterPlot(spikeTimes,stimulusTimes,textConditions(conditionIndices),[-0.05 0 0.2],[],true,0.01,[],{'Condition'});

            set(rasterFig,'Position',[0 0 1920 1080]);
            figfile = sprintf('%s\\multiraster_%s_channel_%s_cluster_%d',stimDir,stimFile,channel,cluster);
            saveas(rasterFig,figfile,'fig');
            saveas(rasterFig,figfile,'png');
        end
        
        nSpikes = zeros(nTrials,nConditions);
        trials = ones(nConditions,1);
        
        for ii = 1:nStimuli
            condition = conditionIndices(ii);
            nSpikes(trials(condition),condition) = sum(spikeTimes > stimulusTimes(ii) & spikeTimes <= stimulusTimes(ii)+0.15);
            trials(condition) = trials(condition)+1;
        end
        
        allNSpikes(:,:,end+1) = nSpikes;
        cells(end+1,:) = [str2double(channel) cluster];
        
        quartiles = prctile(nSpikes,[50 25 75]);
        quartiles = [quartiles(1,:); quartiles(1,:)-quartiles(2,:); quartiles(3,:)-quartiles(1,:)];
        
        if strcmp(options.figs,'yes') 
            figure;
            set(gcf,'Position',[0 0 1200 900]);
            yy = [0 0];
            
            for ii = 1:numel(pausePWs)
                subplot(rows,cols,ii);
                hold on;
                
                legendEntries = {};
                
                for jj = 1:numel(relPWs)
                    indices = conditions(:,2) == relPWs(jj) & conditions(:,3) == pausePWs(ii);
                    
                    if sum(indices) == 0
                        continue;
                    end
                    
                    errorbar( ...
                        conditions(indices,1),                  ...
                        quartiles(1,indices),                   ...
                        quartiles(2,indices),                   ...
                        quartiles(3,indices),                   ...
                        'Color', colours(jj-(jj-1)*(ii==1),:)   ...
                        );
                    
                    legendEntries{end+1} = sprintf('%3.2f%%',relPWs(jj)); %#ok<AGROW>
                end
                
                yy = max(yy,ylim);
                
                if ii > 1
                    legend(legendEntries,'Location','SouthEast');
                    title(sprintf('Pause Time = %d ms',100*pausePWs(ii)));
                else
                    title('No Pause');
                end
            end
            
            for ii = 1:numel(pausePWs)
                subplot(rows,cols,ii);
                
                set(gca,'XScale','log','XTick',[5 10 20 40 80]);
                
                xlim([4 100]);
                ylim(yy);
                
                if ii == (rows-1)*cols+1
                    xlabel('Total Pulse Width/ms');
                    ylabel('# Spikes');
                end
            end
            
            figfile = sprintf('%s\\multispikes_%s_channel_%s_cluster_%d',stimDir,stimFile,channel,cluster);
            saveas(gcf,figfile,'fig');
            saveas(gcf,figfile,'png');
            close(gcf);
        end
    end
    
    forEachChannel(stimFile,[],true,@fn);
    
    save(sprintf('%s\\%s_multi_responses',stimFile,stimDir),'allNSpikes','cells');
end

function uledStimuliToPSRH(doPSRH,stimDir,stimFile,stimuli,ttls,responsiveCells,stimulus,fields,factors,prefix,opts)
    nFields = numel(fields);
    levels = zeros(1,nFields);
    valuess = cell(1,nFields);
    groups = cell(1,nFields-1);

    periodIndex = [];
    paramss = stimuli.(stimulus);
    for hh = 1:nFields
        if ischar(fields{hh}) && strcmp(fields(hh),'T')
            periodIndex = hh;
        end

        if iscell(fields{hh})
            values = zeros(numel(paramss),numel(fields{hh}));

            for gg = 1:numel(fields{hh})
                values(:,gg) = vertcat(paramss.(fields{hh}{gg}));
            end
        else
            values = {paramss.(fields{hh})};
            values = values(:);

            if ~ischar(values{1})
                values = cell2mat(values);
            end
        end

        if hh == 1
            conditions = values;
        else
            groups{hh-1} = values;
        end

        valuess{hh} = unique(values,'rows');
        levels(hh) = size(valuess{hh},1);
    end

    % TODO : mixed array/cell array groups?
    groups = [groups{:}];
    
    cells = zeros(0,3);

    if doPSRH
        options = struct('getopt__version',1,'proportionalspiketimes',false,'rasterfilesuffix',false,'bw',0.1);
        
        uniqueGroups = unique(groups,'rows');

        rasters = cell([0 levels 2]);
        stimulusTimings = cell([levels 2]);
        stimulusLengths = zeros([levels 2]);
    %     subStimulusLengths = zeros(8,2,1,2);
        histograms = cell([0 levels 2]);
        edgess = cell([0 levels 2]);
        repeats = zeros([0 levels 2]);

        meanDiffTimings = 0;
        n = 0;
        for hh = 1:numel(paramss)
            params = paramss(hh);
            subscript = cell(1,nFields+1);

            for gg = 1:nFields
                field = fields{gg};

                if iscell(field)
                    value = zeros(1,numel(field));

                    for ff = 1:numel(field)
                        value(ff) = params.(field{ff});
                    end

                    index = find(ismember(valuess{gg},value,'rows'));
                elseif ischar(params.(field))
                    index = find(strcmp(params.(field),valuess{gg}));
                elseif islogical(params.(field))
                    index = params.(field)+1;
                else
                    index = find(params.(field) == valuess{gg});
                end

                subscript{gg} = index;
            end

            if numel(stimuli.b) == numel(params)
                if stimuli.b(hh).index > numel(ttls)
                    timings = ttls(params.index)+[0 meanDiffTimings/n];
                else
                    timings = ttls([params.index; stimuli.b(hh).index]);
                end
            elseif params.index == numel(ttls)
                timings = ttls(params.index) + [0 params.T/1000];
            else
                timings = ttls(params.index + [0 1]);
            end

            n = n + 1;
            meanDiffTimings = meanDiffTimings + diff(timings);

            subscript{end} = 1;
            index = sub2ind(size(stimulusTimings),subscript{:});

            if isempty(stimulusTimings{index})
                stimulusTimings{index} = timings;
            else
                stimulusTimings{index}(:,1,end+1) = timings;
            end

            stimulusLength = diff(timings);
            stimulusLengths(index) = max(stimulusLengths(index),stimulusLength);
        end

        subStimulusLengths(:,:,1,:) = stimulusLengths;
    end
    
    commandIndices = [stimuli.(stimulus).index]';
    
    doSTA = false;
    
    if isfield(opts,'pixels')
        pixels = opts.pixels;
    
        if ~isempty(pixels) && ~any(isnan(reshape(pixels,numel(pixels),1)))
            pixelsSize = size(pixels);
            lastNonSingletonDimension = find(levels > 1,1,'last');
            nDims = ndims(pixels);
            doSTA = nDims == lastNonSingletonDimension+3 && isequal(pixelsSize(4:lastNonSingletonDimension+3),levels(1:lastNonSingletonDimension));
        end
    end
    
    if ~doPSRH && ~doSTA
        return;
    end
    
    if strcmp(opts.figs,'yes')
        staFig = figure;
        set(staFig,'Visible','off');
    end
    
    if doSTA
        stas = {};
    end
    
    function fn(spikeTimes,~,channel,cluster,~,~)
        if ~isempty(responsiveCells) && ~ismember([str2double(channel) cluster],responsiveCells,'rows')
            return;
        end
        
        tic;
        
        cells = [cells; double(channel) cluster];
        
        if doSTA
            sta = zeros([16 16 pixelsSize(5:end)]);

            for mm = 1:numel(commandIndices)
                nIndices = ndims(pixels)-3;
                indices = cell(1,nIndices);
                
                if iscell(valuess{1})
                    indices{1} = strcmp(conditions{mm},valuess{1});
                else
                    indices{1} = ismember(valuess{1},conditions(mm,:),'rows');
                end
                
                for nn = 2:nIndices
                    indices{nn} = ismember(valuess{2},groups(mm,nn-1));
                end
                
                if numel(periodIndex) == 1 && periodIndex > 1
                    rateMultiplier = 1000/groups(mm,periodIndex-1);
                else
                    rateMultiplier = 1;
                end

                for nn = 1:size(pixels,3)
                    if commandIndices(mm)+nn > numel(ttls)
                        interval = ttls(commandIndices(mm)+nn-1)+[0 paramss(mm).T/1000];
                    else
                        interval = ttls(commandIndices(mm)+nn-[1 0]);
                    end

                    nSpikes = sum(spikeTimes >= interval(1) & spikeTimes < interval(2));
                    rate = nSpikes*rateMultiplier;

                    sta(:,:,indices{2:end}) = sta(:,:,indices{2:end}) + rate*pixels(:,:,nn,indices{:});
                end
            end
            
            stas{end+1} = sta; %#ok<SETNU>
            
            if strcmp(opts.figs,'yes')
                figure(staFig);
                
                figFilePrefix = sprintf('%s\\%s_sta_%s_channel_%s_cluster_%d',stimFile,prefix,stimFile,channel,cluster);
                staSize = size(sta);
                
                for mm = 1:prod(staSize(3:end))
                    clf;
                    surf(sta(:,:,mm));
                    view(2);
                    xlim([1 16]);
                    ylim([1 16]);
                    shading flat;
                    
                    figFile = figFilePrefix;
                    for nn = 1:ndims(sta)-2
                        rowIndex = ind2sub(staSize(3:end),mm);
                        figFile = sprintf('%s_%s_%d',figFile,factors{nn+1},uniqueGroups(rowIndex,nn));
                    end
                        
                    saveas(staFig,figFile,'fig');
                    saveas(staFig,figFile,'png');
                end
            end
        end
        toc;
        
        if ~doPSRH
            return;
        end
        
        [figs,lines,edges,hists,trials] = rasterPlot(spikeTimes,ttls(commandIndices),conditions,[-0.25 0 1.4 1.65],[],true,0.1,groups,factors);
        
        rasters = [rasters; cat(ndims(rasters),reshape(lines,[1 size(lines)]),cell([1 size(lines)]))];
        edgess = [edgess; cat(ndims(edgess),reshape(repmat({edges},levels),[1 levels]),cell([1 levels]))];
        
        dimensionDistributions = cell(1,nFields);
        
        for ii = 1:nFields
            dimensionDistributions{ii} = ones(levels(ii),1);
        end
        
        histSize = size(hists);
        hists = [hists; zeros([1 histSize(2:end)])]; % cpsrh makes histograms with the same number of bins as edges, whereas rasterPlot has one less bin.  Hence pad out the extra bin with zeros (as rasterPlot will never put anything in the last bin)
        histograms = [histograms; cat(ndims(histograms),mat2cell(hists,numel(edges),dimensionDistributions{:}),cell([1 levels]))];
        
        repeats = [repeats; cat(ndims(repeats),reshape(trials,[1 size(trials)]),zeros([1 size(trials)]))];
        
        if ~strcmp(opts.figs,'yes')
            return;
        end
        
        figPrefix = sprintf('%s\\%s_raster_%s_channel_%s_cluster_%d',stimDir,prefix,stimFile,channel,cluster);
        
        for ii = 1:numel(figs)
            set(figs(ii),'Position',[0 0 1600 900]); %,'Visible','on');
            
            figFile = figPrefix;
            for jj = 2:nFields
                figFile = sprintf('%s_%s_%d',figPrefix,factors{jj},uniqueGroups(ii,jj-1));
            end
            
            try
                saveas(figs(ii),figFile,'fig');
                saveas(figs(ii),figFile,'png');
                close(figs(ii));
            catch err
                warning('Unable to save %s raster number %d for channel %s cluster %d in recording %s, skipping...\n',prefix,ii,channel,cluster,stimFile);
                logMatlabError(err);
                continue;
            end
        end
%     (spikeTimes,stimulusTimes,conditions,tmin,tmax,isPSTH,bw,groups,varNames)
    end

    forEachChannel(stimFile,[],true,@fn);
    
    if doPSRH
        save(sprintf('%s\\%s_%s_psrh',stimDir,stimFile,prefix),'options','cells','factors','levels','valuess','rasters','stimulusTimings','stimulusLengths','subStimulusLengths','histograms','edgess','repeats');
    end
    
    if doSTA
        cells = [str2num(char(cells(:,1:2))) cells(:,3)]; %#ok<ST2NM>
        save(sprintf('%s\\%s_%s_sta',stimDir,stimFile,prefix),'cells','stas','factors','conditions','groups');
    end
end

function nothing = movingBars(doPSRH,stimFile,stimDir,stimuli,ttls,spikeRates,lambdas,opts,responsiveCells)
    nothing = NaN;
    
    % TODO : multiple widths?
    width = stimuli.m(1).W;
    pathLength = 16-width;

    pixels = zeros(16,16,pathLength,8);

    [Y,X] = ndgrid(0:15,0:15);
        
    angles = [4; 3; 2; 1; 5; 6; 7; 0];

    for ii = 1:8
        theta = angles(ii)*pi/4;

        for jj = 0:pathLength-1
            w = (X-8)*cos(theta)-(Y-8)*sin(theta)+8;
            pixels(:,:,jj+1,ii) = (w >= jj & w < jj + width);
        end
    end

    opts.pixels = repmat(pixels,[1 1 1 1 2]);

    uledStimuliToPSRH(true,stimDir,stimFile,stimuli,ttls,responsiveCells,'m',{'D' 'T'},{'direction' 'period'},'moving_bar',opts);
end

function nothing = rings(stimFile,stimDir,stimuli,ttls,spikeRates,lambdas,opts,responsiveCells)
    nothing = NaN;
    uledStimuliToPSRH(true,stimDir,stimFile,stimuli,ttls,responsiveCells,'e',{'isContracting' 'T'},{'direction' 'period'},'rings',opts);
end

function nothing = circles(stimFile,stimDir,stimuli,ttls,spikeRates,lambdas,opts,responsiveCells)
    nothing = NaN;
    
    circs = unique([vertcat(stimuli.e.X) vertcat(stimuli.e.Y) vertcat(stimuli.e.R) vertcat(stimuli.e.W)],'rows');
    nCircs = size(circs,1);
    
    pixels = zeros(16,16,1,nCircs);
    [Y,X] = ndgrid(0:15,0:15);
    
    for ii = 1:nCircs
        S = (X-circs(ii,1)).^2 + (Y-circs(ii,2)).^2;
        pixels(:,:,1,ii) = S <= circs(ii,3)^2 & S >= 0;
    end

    opts.pixels = pixels;
    
    uledStimuliToPSRH(true,stimDir,stimFile,stimuli,ttls,responsiveCells,'e',{{'X' 'Y' 'R' 'W'} 'T'},{{'X' 'Y' 'R' 'W'} 'period'},'circle',opts);
end

function nothing = squarePSRH(doPSRH,stimFile,stimDir,stimuli,ttls,spikeRates,lambdas,opts,responsiveCells)
    nothing = NaN;
    
    squares = unique([vertcat(stimuli.r.X) vertcat(stimuli.r.Y)],'rows');
    nSquares = size(squares,1);
    widths = unique(vertcat(stimuli.r.W));
    nWidths = numel(widths);
    
    % TODO : different pulse widths
    pixels = zeros(16,16,1,nSquares,nWidths);
    
    for hh = 1:nWidths
        w = widths(hh);

        for ii = 1:nSquares
            if squares(ii,1)+w > 16 || squares(ii,2)+w > 16
                continue;
            end
            
            pixels(squares(ii,2)+(1:w),squares(ii,1)+(1:w),1,ii,hh) = 1;
        end
    end

    opts.pixels = pixels;
    
    uledStimuliToPSRH(doPSRH,stimDir,stimFile,stimuli,ttls,responsiveCells,'r',{{'X' 'Y'} 'W' 'T'},{{'X' 'Y'} 'width' 'period'},'square',opts);
end

function nothing = shapes(stimFile,stimDir,stimuli,ttls,spikeRates,lambdas,opts,responsiveCells)
    nothing = NaN;
    
    nStimuli = numel(stimuli.i);
    
    commandIndices = [stimuli.i.index]';
    T = [stimuli.i.T]';
    F = {stimuli.i.F}';
    
    pws = unique(T);
    maxPW = max(pws);
    responseWindow = 2*maxPW/1000;
   
    shapes = unique(F);
    
    I = zeros(size(F));
    
    for ii = 1:nStimuli
        I(ii) = find(strcmp(F{ii},shapes));
    end
    
    nPWs = numel(pws);
    nShapes = numel(shapes);
    nConditions = nPWs*nShapes;
    nTrials = nStimuli/nConditions;
    
    nSpikes = zeros(nTrials,nShapes,nPWs,0);
    nCells = 0;
    
    cells = zeros(0,2);
    
    function fn(spikeTimes,~,channel,cluster,~,~)
        if ~isempty(responsiveCells) && ~ismember([str2double(channel) cluster],responsiveCells,'rows')
            return;
        end
        
        nCells = nCells + 1;
        
        cells(end+1,:) = [str2double(channel) cluster]; %#ok<SETNU>
        
        [figs,lines] = rasterPlot(spikeTimes,ttls(commandIndices),I,[-0.05 0 responseWindow+[0 0.05]],[],true,0.01,T,{'Image' 'Pulse Width'});
        
        for jj = 1:nShapes
            for kk = 1:nPWs
                raster = lines{jj,kk};
                
                for ll = 1:nTrials
                    spikes = raster{ll};
                    nSpikes(ll,jj,kk,nCells) = sum(spikes > 0 & spikes <= responseWindow);
                end
            end
        end
        
        if ~strcmp(opts.figs,'yes')
            return;
        end
        
        figPrefix = sprintf('%s\\shapes_raster_%s_channel_%s_cluster_%d',stimDir,stimFile,channel,cluster);
        
        for jj = 1:numel(figs)
            set(figs(jj),'Position',[0 0 1600 900]);
            figFile = sprintf('%s_pw_%d',figPrefix,pws(jj));
            saveas(figs(jj),figFile,'fig');
            saveas(figs(jj),figFile,'png');
        end
    end

    forEachChannel(stimFile,[],true,@fn);
    
    maxNSpikes = max(max(max(max(nSpikes))));
    
    % p(s|r) = p(r|s)*p(s)/p(r) - the way we've done the randomisation
    % ensures p(s) is constant for all s and p(r) is fixed for a given
    % response, so to maximise p(s|r) over s we need only maximise p(r|s)
    respStimCPDF = zeros(maxNSpikes+1,nCells,nShapes,nPWs);
    
    randomIndices = randperm(nTrials);
    trainIndices = randomIndices(1:ceil(nTrials/2));
    
    for ii = 1:nShapes
        for hh = 1:nPWs
            for gg = 1:nCells
                respStimCPDF(:,gg,ii,hh) = hist(nSpikes(trainIndices,ii,hh,gg),0:maxNSpikes);
            end
        end
    end
       
    respStimCPDF = respStimCPDF./repmat(sum(respStimCPDF),[maxNSpikes+1 1 1 1]);
    respStimCPDF = safelog(respStimCPDF);
    
    performance = zeros(nPWs,1);
    testIndices = randomIndices(1+ceil(nTrials/2):end);
    
    for ii = 1:nPWs
        for hh = 1:nShapes
            for gg = 1:numel(testIndices)
                index = testIndices(gg);
                counts = squeeze(nSpikes(index,hh,ii,:));
                posterior = zeros(nShapes,1);
            
                for ff = 1:nCells
                    posterior = posterior + squeeze(respStimCPDF(counts(ff)+1,ff,:,ii));
                end
                
                [~,maxS] = max(posterior);
                
                if maxS == hh
                    performance(ii) = performance(ii) + 1;
                end
            end
        end
    end    
    
    performance = performance/(nShapes*nTrials/2); %#ok<NASGU>
    
    save(sprintf('%s\\%s_shape_decoder',stimDir,stimFile),'cells','nSpikes','respStimCPDF','trainIndices','testIndices','performance','pws');
end