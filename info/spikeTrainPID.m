function pid = spikeTrainPID(s,R,labels)
    tic;
    [N,nSources] = size(R);
    
    assert(N == numel(s),'There must be exactly one response observation for each stimulus presentation.');
    
    stimuli = unique(s);
    nStimuli = numel(unique(s));
    
    kStimuli = zeros(nStimuli,1);
    
    for ii = 1:nStimuli
        kStimuli(ii) = sum(s == stimuli(ii));
    end
    
    pStimuli = kStimuli/N;
    
    individualInformation = zeros(nSources,1);
    specificInformation = zeros(nStimuli,nSources); % planes : non-info & info singletons, hyperplanes : raw, max_embed_dim 2, max_embed_dim nMax
    
    nSpikes = cellfun(@numel,R);
    nMax = max(nSpikes(:));
    E = zeros(N,nMax+1,nSources);
    countss = zeros(N,nSources);
    
    toc;
    
    for ii = 1:nSources
        tic;
        r = R(:,ii);
        
        ns = unique(nSpikes(:,ii));
        n = numel(ns);
        
        Z = struct('M',int32(8),'N',int32(1));
        Z.sites = struct('label',{labels(ii)},'recording_tag',{{'episodic'}},'time_scale',1,'time_resolution',1/7022,'si_unit','none','si_prefix',1);
        Z.categories = struct('label',{{'1'};{'2'};{'3'};{'4'};{'5'};{'6'};{'7'};{'8'}},'P',int32(125));
        for jj = 1:8
            stimulusIndices = s == jj;
            spikeTrains = r(stimulusIndices);
            Z.categories(jj).trials = struct('start_time',0,'end_time',0.5,'Q',cellfun(@(x) int32(numel(x)),spikeTrains,'UniformOutput',false),'list',spikeTrains);
        end
        
        [data,counts,categories] = binlessopen(Z);
        countss(:,ii) = counts;
        
        warped = binlesswarp(data);
        embedded = binlessembed(warped,struct('max_embed_dim',nMax));
        E(:,:,ii) = embedded;

        X = cell(nStimuli,n);
        Y = cell(1,n);
        
        for jj = 1:n
            hasNSpikes = nSpikes(:,ii) == ns(jj);
                
            Y{jj} = embedded(hasNSpikes,2:end);
                
            for kk = 1:nStimuli
                X{kk,jj} = embedded(hasNSpikes & s == stimuli(kk),2:end);
            end
        end
        
        degenerate = ns == 0;
                    
        for jj = 1:nStimuli
            specificInformation(jj,ii) = mixedKLDivergence(X(jj,:),Y,degenerate,1,'both',true);
        end
            
        X = cell(1,n);
        
        for jj = 1:n
            hasNSpikes = nSpikes(:,ii) == ns(jj);
            X{jj} = s(hasNSpikes);
        end
        
        [~,~,~,I] = binlessinfo(embedded,counts,categories,Z.M,struct('stratification_strategy',1,'singleton_strategy',1));
        individualInformation(ii) = I.value;
            
        fprintf('Calculated KLD and MI for neuron %d of %d in %f seconds\n',ii,nSources,toc);
    end
    
    nPairs = nSources*(nSources-1)/2;
    
    pairs = zeros(nPairs,2);
    ensembleInformation = zeros(nPairs,1);
    
    n = 0;
    for ii = 1:nSources-1
        for jj = (ii+1):nSources
            tic;
            n = n + 1;
            pair = [ii jj];
            pairs(n,:) = pair;
            embedded = E(:,:,pair);
            counts = countss(:,pair);
            [~,~,~,I] = binlessinfom(embedded,counts,categories,Z.M,struct('stratification_strategy',1,'singleton_strategy',1));
            ensembleInformation(n) = I;
            fprintf('Calculated MI for pair %d of %d in %f seconds\n',n,nPairs,toc);
        end
    end
    
    tic;
    redundancy = (sum(repmat(pStimuli,1,n).*min(cat(3,specificInformation(:,pairs(:,1)),specificInformation(:,pairs(:,2))),[],3)))';
    
    uniqueInformation = [individualInformation(pairs(:,1)) individualInformation(pairs(:,2))]-repmat(redundancy,1,2);
    
    synergy = ensembleInformation-individualInformation(pairs(:,1))-individualInformation(pairs(:,2))+redundancy;
    
    pid = [redundancy uniqueInformation synergy ensembleInformation];
    toc;
end