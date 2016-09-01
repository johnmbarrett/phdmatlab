function pids = biasCorrectedPID(x,Y,varargin)
    options = getopt('savefile=NaN',varargin{:});
    
    sizeX = size(x);
    [s,~,x] = unique(x);
    ns = numel(s);
    nx = sizeX(1);

    if iscell(Y)
        d = ndims(Y{1});
        
        if size(Y,2) == 1
            Y = permute(cat(d+1,Y{:}),[d+1 1:d]);
            isShuffle = false;
        elseif size(Y,2) == 2
            Y = permute(cat(d+2,cat(d+1,Y{:,1}),cat(d+1,Y{:,2})),[d+1 1:d d+2]);
            isShuffle = true;
        else
            error('I don''t know what you expect me to do here.');
        end
    end
    
    sizeY = size(Y);
    assert(nx == sizeY(2));
    n = sizeY(1);
    N = n*(n-1)/2;
    M = prod(sizeY(3:end-isShuffle));
    assert(prod(sizeX(2:end)) == M);
    
    pairs = zeros(N,2);
    indices = zeros(n);
    
    nn = 0;
    for ii = 1:n-1
        for jj = (ii+1):n
            nn = nn + 1;
            indices(ii,jj) = nn;
            pairs(nn,:) = [ii jj];
        end
    end
    
    kld = zeros(n,M,ns);
    mi1 = zeros(n,M);
    mi2 = zeros(N,M);

    x = reshape(x,nx,M);
    Y = reshape(Y,[n nx M 1+isShuffle]);
    ps = zeros(1,M,ns);
    
    for ii = 1:M
        ps(1,ii,:) = accumarray(x(:,ii),1)/nx;
    end

    for hh = 1:M
        for ii = 1:n
            tic;
            if isShuffle
                [~,~,y] = unique([Y(ii,:,hh,1) Y(ii,:,hh,2)]');
                Y(ii,:,hh,1) = y(1:nx);
                Y(ii,:,hh,2) = y(nx+1:2*nx);
                y = y(1:nx);
            else
                [~,~,y] = unique(Y(ii,:,hh,1)');
                Y(ii,:,hh,1) = y;
            end
            
            for jj = 1:ns
                kld(ii,hh,jj) = biasCorrect(y(x(:,hh) == jj),y,'method','function','mifunction',@(x,y,varargin) discreteKLD(x,y));
            end
            
            mi1(ii,hh) = biasCorrect(x(:,hh),y,'method','discrete');
            toc;
        end
        
        for ii = 1:n-1
            for jj = (ii+1):n
                tic;
                if isShuffle
                    yy = [Y(ii,:,hh,1)' Y(jj,:,hh,2)'];
                else
                    yy = Y([ii jj],:,hh,1)';
                end
                
                mi2(indices(ii,jj),hh) = biasCorrect(x(:,hh),yy,'method','discrete');
                toc;
            end
        end
    end
    
    if ischar(options.savefile)
        save(options.savefile,'-v7.3','kld','mi1','mi2','ps','pairs','indices','x','Y','sizeX','sizeY');
    end

    pids = zeros(N,M,5);
    
    pids(:,:,1) = sum(repmat(ps,N,1).*min(cat(4,kld(pairs(:,1),:,:),kld(pairs(:,2),:,:)),[],4),3);
    pids(:,:,2) = mi1(pairs(:,1),:)-pids(:,:,1);
    pids(:,:,3) = mi1(pairs(:,2),:)-pids(:,:,1);
    pids(:,:,5) = mi2;
    pids(:,:,4) = pids(:,:,5) - sum(pids(:,:,1:4),3);
    
    pids = reshape(permute(pids,[1 3 2]),[N 5 sizeY(3:end-isShuffle)]);
    
    if ischar(options.savefile)
        save(options.savefile,'-append','pids');
    end
end