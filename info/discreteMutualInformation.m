function I = discreteMutualInformation(x,y,varargin)
    [pxy,nx,ny] = ejpmf(x,y);
    
    px = epmf(x);
    py = epmf(y);
    
    pxpy = repmat(px,1,ny).*repmat(py',nx,1);
    
    I = sum(sum(plogq(pxy,pxy./pxpy)));
    
    options = getopt('base=2',varargin{:});
    
    I = I/log(options.base);
end