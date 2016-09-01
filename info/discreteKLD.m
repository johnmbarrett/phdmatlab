function D = discreteKLD(x,y,varargin)
    % TODO : multiple dimensions
    assert(isvector(x) & isvector(y));
    
    x = x(:);
    y = y(:);
    
    options = getopt('base=2 edges=[] isPMF=false',varargin{:});
    edges = options.edges;
    isPMF = options.ispmf;
    
    if isempty(isPMF) || ~all(logical(isPMF(:)))
        if ~isempty(edges) && isnumeric(edges) && all(isfinite(edges(:)))  % assume samples from discrete distributions
            if isscalar(edges)
                edges = linspace(min([x;y]),max([x;y]),edges);
            else
                if ~isvector(edges);
                    edges = edges(:);
                end

                if ~issorted(edges)
                    edges = sort(edges);
                end
            end

            x = histc(x,edges);
            x = x/sum(x);

            y = histc(y,edges);
            y = y/sum(y);
        else
            [~,~,idx] = unique([x;y]);
            
            nx = numel(x);
            ix = idx(1:nx);
            x = accumarray(ix,1)/nx;
            
            iy = idx((nx+1):end);
            y = accumarray(iy,1)/numel(y);
            
            ux = unique(ix);
            uy = unique(iy);
            
            x(setdiff(uy,ux),1) = 0;
            y(setdiff(ux,uy),1) = 0;
        end
    end
    
    if any(y == 0 & x ~= 0)
        D = Inf;
        return;
    end
    
    D = x.*(log(x)-log(y));
    
    D(x == 0) = 0;
    
    D = sum(D);
    
    D = D/log(options.base);
end