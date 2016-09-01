function D = kdeKLD(p,q,varargin)
    options = getopt('n=2^12 base=2 minx=NaN maxx=NaN',varargin{:});
    
    if ~isnan(options.minx)
        minX = options.minx;
    else
        minX = min(min(p(:)));
    end
    
    if ~isnan(options.maxx)
        maxX = options.maxx;
    else
        maxX = max(max(q(:)));
    end
    
    [~,pc,x] = kde(p,options.n,minX,maxX);
    [~,qc,y] = kde(q,options.n,minX,maxX);
    
    pc = pc/trapz(x,pc);
    qc = qc/trapz(y,qc);
    
    if ~isequal(x,y)
        qc = interp1(y,qc,x);
        
        if any(isnan(qc(:)))
            D = Inf;
            return;
        end
    end
    
    D = plogq(pc,pc./qc);
    
    if any(isinf(D(:)))
        D = Inf;
        return;
    end
    
    D(isnan(D)) = 0;
    
    D = trapz(x,D);
    
    D = D/options.base;
end