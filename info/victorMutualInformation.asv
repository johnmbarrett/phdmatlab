function [I,Icount,Itiming] = victorMutualInformation(S,R,degenerate,includeSingletons,useVictorEmbedding)
    if nargin < 5
        useVictorEmbedding = false;
    end
    
    if ~iscell(S)
        S = {S};
    end
    
    if ~iscell(R)
        R = {R};
    end
    
    Nn = cellfun(@numel,S);
    N = sum(Nn);
    
    assert(isequal(Nn,cellfun(@(A) size(A,1),R)));
    
    nMax = numel(S);
    
    degenerate = logical(degenerate) | Nn == 1;
    assert(nMax == numel(degenerate));
    
    Itiming = zeros(nMax,1);
    
    indices = 1:nMax;
    
    for ii = indices(~degenerate)
        Itiming(ii) = victorContinuousMutualInformation(S{ii},R{ii},includeSingletons,useVictorEmbedding);
    end
    
    Icount = -sum(plogp(Nn));
    
    Iak = cell(nMax,1);
    Jak = cell(nMax,1);
    Nak = cell(nMax,1);
    
    sMax = -Inf;
    allStimuli = [];
    
    for ii = indices
        s = S{ii};
        stimuli = unique(s);
        nStimuli = numel(stimuli);
        
        if nStimuli == 0
            continue;
        end
           
        allStimuli = union(stimuli,allStimuli);
        
        sMax = max(sMax,max(stimuli));
        
        Jak{ii} = stimuli(:);
        
        Iak{ii} = ii*ones(nStimuli,1);
        
        Nak{ii} = zeros(nStimuli,1);
        
        for jj = 1:nStimuli
            Nak{ii}(jj) = sum(s == stimuli(jj));
        end
        
        Icount = Icount + sum(plogp(Nak{ii}));
    end
    
    Icount = Icount/N;
    
    qk = full(sum(sparse(vertcat(Iak{:}),vertcat(Jak{:}),vertcat(Nak{:}),nMax,sMax),1))/N;
    
    Icount = Icount - sum(plogp(qk)) - (numel(allStimuli)-1)*(nMax-1)/(2*N*log(2));
    
    Itiming = sum(Nn(:).*Itiming/N);
    I = Icount + Itiming;
end

function I = victorContinuousMutualInformation(s,R,includeSingletons,legendreCoeffs)
    assert(isvector(s),'Stimulus variable must be a vector');
    assert(ismatrix(R),'Response variable must be a matrix');
    
    [N,d] = size(R);
    
    s = s(:);
    
    assert(numel(s) == N,'There must be one response observation for every stimulus');
    
    if ~isscalar(legendreCoeffs) && istriu(legendreCoeffs) && size(legendreCoeffs,1) == size(legendreCoeffs,2)
        J = tiedrank(R(:));
        R = 2*((reshape(J,size(R))-0.5)/(N*d))-1;

        r = min(d,size(legendreCoeffs,1));
        
        P = zeros(N,r);
        
        for ii = 1:r
            % TODO : quicker to precompute the legendre coefficients and
            % use polyval or don't compute coeffs and use the recurrence?
            P(:,ii) = sqrt(2*ii+1)*sum(polyval(legendreCoeffs(1:ii+1,ii+1),R),2);
        end
        
        R = P;
        d = r;
        clear P;
    end
    
    [s,sortIndices] = sort(s);
    R = reshape(R(sortIndices,:),[N 1 d]);
    
    rho = sqrt(sum((repmat(R,[1 N 1])-repmat(permute(R,[2 1 3]),[N 1 1])).^2,3));
    
    [rhoSorted,sortIndices] = sort(rho,2);
    lambda = rhoSorted(:,2);
    
    C = lambda ~= 0;
    
    [stimuli,~,sidx] = unique(s);
    
    S = numel(stimuli);
    
    Na = zeros(S,1);
    
    for ii = 1:S
        Na(ii) = sum(s == stimuli(ii));
    end
    
    singletons = false(S,1);
    
    % There's a bizarre edge case where all but one response to at least
    % one stimulus is part of a zero-distance set, so by Na that response
    % is not a singleton, but it has no counterpart when you calculate
    % per-stimulus nearest-neighbour distances for the set of distinct
    % responses and hence *is* a singleton in that sense
    for ii = 1:S
        singletons(ii) = sum(s(C) == stimuli(ii)) == 1;
    end
    
    singletonStimuli = stimuli(singletons);
    
    C = C & ~ismember(s,singletonStimuli);
    Nc = sum(C);
    Is = cell(S,1);
    Ns = zeros(S,1);

    if Nc == 0
        Icontinuous = 0;
    else
        for ii = 1:S
            idx = find(C & s == stimuli(ii));
            Is{ii} = idx;
            Ns(ii) = numel(idx);
        end

        assert(all(Na >= Ns));

        cNs = [0; cumsum(Ns)];

        assert(cNs(end) == Nc);
        assert(~any(diff(vertcat(Is{:})) <= 0));

        nu = Inf(Nc,max(Ns));

        for ii = 1:numel(Ns)
            nu(cNs(ii)+1:cNs(ii+1),1:Ns(ii)) = rho(Is{ii},Is{ii});
        end

        nu = sort(nu,2);
        mu = nu(:,2);

        assert(all(mu > 0));

        Icontinuous = (d/N)*sum(log(lambda(C)./mu))-sum((Ns/Nc).*log((Ns-1)/(Nc-1)));
    
        if all(C)
            I = Icontinuous;
            return;
        end
    end
        
    Z = ~C;
    
    zidx = find(Z);
    Nz = zeros(numel(zidx),1);
    Nzs = zeros(numel(zidx),S);
    Zn = zeros(N,1);
    
    n = 0;
    for ii = 1:numel(zidx)
        row = zidx(ii);
        
        if Zn(row) > 0
            continue;
        end
        
        % TODO : is using find(...,1) quicker here?
        cols = sortIndices(row,rhoSorted(row,:) == 0);
        n = n + 1;
        
        Zn(cols) = n;
        
        nCols = numel(cols);
        Nz(n) = nCols;
        
        for jj = 1:nCols
            Nzs(n,sidx(cols(jj))) = Nzs(n,sidx(cols(jj))) + 1;
        end
    end
    
    assert(sum(Nz) == numel(zidx));
    
    Nz = [Nz(1:n); Nc];
    Nzs = [Nzs(1:n,:); Ns'];
    
    assert(sum(Nz) == N);
    assert(all(sum(Nzs)' == Na));
    
    if ~includeSingletons
        nonSingletonSets = Nz ~= 1;
        Nz = Nz(nonSingletonSets);
        Nzs = Nzs(nonSingletonSets,:);
    end
    
    Ipartition = (sum(sum(plogp(Nzs)))-sum(plogp(Nz)))/N - sum(plogp(Na/N)) - (S-1)*n/(2*N*log(2));
    
    I = Ipartition + (Nc/N)*Icontinuous;
end