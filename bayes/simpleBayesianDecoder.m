function [fractionCorrect,perfVTrials,perfVUnits] = simpleBayesianDecoder(s,R,subTrials,subUnits,nReps,model,lb,ub,trainPrior)
    if nargin < 9
        trainPrior = false;
    end
    
    assert(numel(s) == size(R,1),'There must be the same number of observations for stimulus and response.');
    
    s = s(:);
    
    valid = ~isnan(s);
    
    if ~iscell(R)
        valid = valid & ~any(isnan(R),2);
    end
    
    
    [u,~,s] = unique(s(valid));
    R = R(valid,:);
    
    N = numel(s);
    m = numel(u);
    k = size(R,2);
    
    if nargin < 7
        if iscell(R)
            lb = min(cellfun(@(r) tertiaryop(isempty(r),Inf,min(r)),R(:)));
        else
            lb = min(R(:));
        end
    end
        
    if nargin < 8
        if iscell(R)
            ub = max(cellfun(@(r) tertiaryop(isempty(r),-Inf,max(r)),R(:)));
        else
            ub = max(R(:));
        end
    end
    
    if nargin < 6
        model = 'discrete';
    end
        
    initDecoder = @(varargin) struct([]);
    
    makeSubDecoderFun = @(S,idx,k) structfun(@(v) makeSubDecoder(v,idx,k),S,'UniformOutput',false);
    
    % this should really be OO at this point
    if isstruct(model)
        assert(isfield(model,'trainDecoder') && isa(model.trainDecoder,'function_handle') && isfield(model,'likelihoodFun') && isa(model.likelihoodFun,'function_handle'),'You must provide a training function and a likelihood function');
        trainDecoder = model.trainDecoder;
        likelihoodFun = model.likelihoodFun;
        
        if isfield(model,'initDecoder') && isa(model.initDecoder,'function_handle')
            initDecoder = model.initDecoder;
        end
        
        if isfield(model,'makeSubDecoder') && isa(model.makeSubDecoder,'function_handle')
            makeSubDecoderFun = model.makeSubDecoder;
        end
    elseif strcmp(model,'kde')
        trainDecoder = @trainKDEDecoder;
        likelihoodFun = @getKDELikelihood;
        initDecoder = @initKDEDecoder;
    elseif strcmp(model,'discrete')
        trainDecoder = @trainDiscreteDecoder;
        likelihoodFun = @getDiscreteLikelihood;
        initDecoder = @initDiscreteDecoder;
    elseif strcmp(model,'mixedBernoulliGaussian')
        trainDecoder = @trainMBGDecoder;
        likelihoodFun = @getMBGLikelihood;
    elseif strcmp(model,'mixedBernoulliKDE')
        trainDecoder = @trainMBKDecoder;
        likelihoodFun = @getMBKLikelihood;
        initDecoder = @initKDEDecoder;
    else
        error('Unknown model %s\n',model);
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
    
    decoderParams = initDecoder(s,R,lb,ub);
    
    nsu = numel(subUnits);
    nst = numel(subTrials);
    
    perfVTrials = zeros(nReps,nst);
    ns = N/m;
    
    if ~trainPrior
        prior = constructPrior(1:N,s,m,N);
    end
    
    for ii = 1:nst
        nt = subTrials(ii);
        
        for jj = 1:nReps
            idx = randperm(ns,nt);
            % TODO : this only works if equal numbers of stimuli
            trainIndices = repmat(idx(:),m,1) + ns*kron((0:m-1)',ones(nt,1));
    
            if trainPrior
                prior = constructPrior(trainIndices,s,m,N);
            end
            
            decoder = trainDecoder(trainIndices,s,R,decoderParams);
            
            testIndices = setdiff(1:N,trainIndices);
            nTests = numel(testIndices);
            
            for kk = 1:nTests
                perfVTrials(jj,ii) = perfVTrials(jj,ii) + testResponse(testIndices(kk),s,R,decoder,prior,likelihoodFun)/nTests;
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
    
        if trainPrior
            prior = constructPrior(trainIndices,s,m,N);
        end
        
%         if exist('newBW','var');
%             oldBW = newBW;
%         end
        
        decoder = trainDecoder(trainIndices,s,R,decoderParams);
        
%         if isfield(decoder,'bw')
%             newBW = decoder.bw;
%             
%             if exist('oldBW','var')
%                 disp(corr(oldBW(:),newBW(:)));
%                 disp(mean(abs(oldBW(:)-newBW(:))));
%             end
%         end
        
        correct(ii) = testResponse(ii,s,R,decoder,prior,likelihoodFun);
        
        for jj = 1:nReps
            for kk = 1:nsu
                idx = unitIndices{jj,kk};
                
                subDecoder = makeSubDecoderFun(decoder,idx,k);
%                 if model
%                     subDecoder = struct('ps',decoder.ps,'prk',decoder.prk(:,idx),'prks',decoder.prks(:,idx,:),'x',decoder.x);
% %                     % this is bad code but it's more efficient for the
% %                     % non-KDE case
% %                     subDecoder = trainDecoder(trainIndices,s,R(:,idx));
%                 else
%                     subDecoder = struct('ps',decoder.ps,'prk',decoder.prk(idx),'prks',decoder.prks(idx));
%                 end
                
                perfVUnits(jj,kk) = perfVUnits(jj,kk) + testResponse(ii,s,R(:,idx),subDecoder,prior,likelihoodFun)/N;
            end
        end
    end
    
    fractionCorrect = sum(correct)/N;
end

function prior = constructPrior(trainIndices,s,m,N)
    prior = zeros(1,m);

    for ii = 1:m
        prior(ii) = sum(s(trainIndices) == ii)/N;
    end
end

function correct = testResponse(testIndex,s,R,model,prior,logLikelihoodFun)
    response = R(testIndex,:);
    logLikelihood = logLikelihoodFun(response,model);
    
    % TODO : something more mathematically rigorous than this
    logLikelihood(isinf(logLikelihood)) = log(eps);
    
    logPosterior = log(repmat(prior,size(R,2),1)) + logLikelihood; % - log(repmat(model,1,size(prior,2)));
    logPosterior = mean(logPosterior,1);

    decodedStimuli = find(logPosterior == max(logPosterior));
    targetStimulus = s(testIndex);

    if numel(decodedStimuli) == 1
        correct = decodedStimuli == targetStimulus;
    elseif ismember(targetStimulus,decodedStimuli)
        correct = 1/numel(decodedStimuli);
    else
        correct = 0;
    end
end