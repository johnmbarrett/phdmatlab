function pids = pairwisePID(x,Y)
    [s,~,x] = unique(x);
    ns = numel(s);
    nx = numel(x);

    if iscell(Y)
        d = ndims(Y{1});
        Y = permute(cat(d+1,Y{:}),[d+1 1:d]);
    end
    
    sizeY = size(Y);
    n = sizeY(1);
    N = n*(n-1)/2;
    M = prod(sizeY(3:end));
    
    pids = zeros(N,4,M);

    [sources,sets,trans] = PIDLattice(2);

    for hh = 1:M
        oo = 1;
        tic;
        for ii = 1:n-1
            for jj = (ii+1):n
                [r1,~,y1] = unique(Y(ii,:,hh));
                [r2,~,y2] = unique(Y(jj,:,hh));
                
                counts = zeros(ns,numel(r1),numel(r2));
                
                for kk = 1:nx
                    ll = x(kk);
                    mm = y1(kk);
                    nn = y2(kk);
                    counts(ll,mm,nn) = counts(ll,mm,nn) + 1;
                end
                
                pids(oo,:,hh) = PID(counts,sources,sets,trans)';
                oo = oo + 1;
            end
        end
        toc;
    end

    pids(:,5,:) = sum(pids,2);
    pids = reshape(pids,[N 5 sizeY(3:end)]);
end