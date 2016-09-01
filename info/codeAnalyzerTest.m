function ret = codeAnalyzerTest(useKDE)
    assert(numel(s) == size(R,1),'There must be the same number of observations for stimulus and response.');
    
    s = s(:);
    
    valid = ~isnan(s) & ~any(isnan(R),2);
    [u,~,s] = unique(s(valid));
    R = R(valid,:);
    
    N = numel(s);
    m = numel(u);
    k = size(R,2);
    
    if nargin >= 1 && useKDE
        if nargin < 7
            lb = min(R(:));
        else
            assert(all(R(:) >= lb),'Empirical response probability must be zero below the support of the response distribution');
        end
        
        if nargin < 8
            ub = max(R(:));
        else
            assert(all(R(:) <= ub),'Empirical response probability must be zero above the support of the response distribution');
        end
        
        trainDecoder = @(trainIndices,s,R) trainKDEDecoder(trainIndices,s,R,lb,ub);
        testResponse = @testKDEResponse;

        return
    else
        useKDE = false;
        vk = cell(k,1);
        nk = zeros(k,1);

        for ii = 1:k
            [v,~,r] = unique(R(:,ii));
            vk{ii} = v;
            nk(ii) = numel(v);
            R(:,ii) = r;
        end
        
        trainDecoder = @(trainIndices,s,R) trainDiscreteDecoder(trainIndices,s,R,vk,nk,m);
        testResponse = @testDiscreteResponse;
    end

    if nargin < 5
        nReps = 0;
    end
    
    if nargin < 4
        subUnits = [];
    end
    
    if nargin < 3
        subTrials = [];
    end
    
    nsu = numel(subUnits);
    nst = numel(subTrials);
    
    perfVTrials = zeros(nReps,nst);
    ns = N/m;
    
    for ii = 1:nst
        nt = subTrials(ii);
        
        for jj = 1:nReps
            idx = randperm(ns,nt);
            % TODO : this only works if equal numbers of stimuli
            trainIndices = repmat(idx(:),m,1) + ns*kron((0:m-1)',ones(nt,1));
            
            [ps,prk,prks] = trainDecoder(trainIndices,s,R);
            
            testIndices = setdiff(1:N,trainIndices);
            nTests = numel(testIndices);
            
            for kk = 1:nTests
                perfVTrials(jj,ii) = perfVTrials(jj,ii) + testResponse(testIndices(kk),s,R,ps,prk,prks)/nTests;
            end
        end
    end
    
    unitIndices = cell(nReps,nsu);
    
    for ii = 1:nsu
        nu = subUnits(ii);
        
        for jj = 1:nReps
            unitIndices{jj,ii} = randperm(k,nu);
        end
    end
    
    perfVUnits = zeros(nReps,nsu);
    
    correct = zeros(N,1);
    
    for ii = 1:N
        trainIndices = setdiff(1:N,ii);
        
        decoder = trainDecoder(trainIndices,s,R);
        
        correct(ii) = testResponse(ii,s,R,decoder);
        
        for jj = 1:nReps
            for kk = 1:nsu
                idx = unitIndices{jj,kk};
                
                if useKDE
                    subDecoder = decoder;
                else
                    subDecoder = struct('ps',decoder.ps,'prk',decoder.prk(idx),'prks',decoder.prks(idx));
                end
                
                perfVUnits(jj,kk) = perfVUnits(jj,kk) + testResponse(ii,s,R(:,idx),ps,prk(idx),prks(idx))/N;
            end
        end
    end
    
    ret = sum(correct)/N;
end

function nb = trainKDEDecoder(trainIndices,s,R,lb,ub)
    nb = NaiveBayes.fit(R(trainIndices,:),s(trainIndices),'Distribution','kernel','KSSupport',[lb ub]);
end

% function [ps,prk,prks] = trainDiscreteDecoder(trainIndices,s,R,vk,nk,m,N)
function pdfs = trainDiscreteDecoder(trainIndices,s,R,vk,nk,m)
    N = numel(trainIndices);
    ps = zeros(1,m);
    prk = cellfun(@(v) zeros(numel(v),1),vk,'UniformOutput',false);
    prks = cellfun(@(v) zeros(numel(v),m),vk,'UniformOutput',false);
    
    for ii = 1:m
        sameStimulus = s(trainIndices) == ii;
        ns = sum(sameStimulus);
        ps(ii) = ns/N;

        for jj = 1:numel(nk)
            n = nk(jj);

            for kk = 1:n
                nrs = sum(sameStimulus & R(trainIndices,jj) == kk);
                prk{jj}(kk) = prk{jj}(kk) + nrs/N;
                prks{jj}(kk,ii) = nrs/ns;
            end
        end
    end
    
    pdfs.ps = ps;
    pdfs.prk = prk;
    pdfs.prks = prks;
end

function correct = testKDEResponse(testIndex,s,R,nb)
    t = predict(nb,R(testIndex,:));
    correct(ii) = t == s(testIndex);
end

function correct = testDiscreteResponse(testIndex,s,R,pdfs)
    ps = pdfs.ps;
    prk = pdfs.prk;
    prks = pdfs.prks;
    

    m = numel(ps);
    k = numel(prk);
    psrk = zeros(m,k);
    testedResponse = R(testIndex,:);
        
    for ii = 1:k
        psrk(:,ii) = (prks{ii}(testedResponse(ii),:)*prk{ii}(testedResponse(ii))./ps)';
    end

    % TODO : something more mathematically rigorous than this
    psrk(psrk == 0) = eps;

    psR = sum(log(psrk),2);

    [~,decodedStimuli] = max(psR);
    targetStimulus = s(testIndex);

    if numel(decodedStimuli) == 1
        correct = decodedStimuli == targetStimulus;
    elseif ismember(targetStimulus,decodedStimuli)
        correct = 1/numel(decodedStimuli);
    else
        correct = 0;
    end
end