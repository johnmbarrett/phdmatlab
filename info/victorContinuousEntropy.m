function H = victorContinuousEntropy(x)
    [N,d] = size(x);
    
    x = reshape(x,[N 1 d]);
    
    lambda = sort(sum((repmat(x,1,N,1)-repmat(permute(x,[2 1 3]),N,1,1)).^2,3),2);
    lambda = lambda(:,2);
    
    S = d*pi.^(d/2)/gamma(d/2+1);
    
    H = (d/N).*sum(log(lambda))+log(S*(N-1)/d)+0.5772156649015328606065120900824;
end