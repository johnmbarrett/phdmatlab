function H = continuousEntropy(x,k)
    if nargin < 2
        k = 1;
    end
    
    [N,d] = size(x);
    
    assert(k < N,'kNN parameter must be stricly less than the number of samples.');
    
    x = reshape(x,[N 1 d]);
    
    rho = sort(sum((repmat(x,1,N,1)-repmat(permute(x,[2 1 3]),N,1,1)).^2,3),2);
    
    V = pi.^(d/2)/gamma(d/2+1);
    
    H = (d/N).*sum(log(rho(:,k+1)))+log(V*(N-1))-psi(k)+log(N)-3;
end