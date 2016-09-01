function Y = fft2reflect(X,replaceinf)
    M = size(X,1);
    N = size(X,2);
    Mdash = 2*M-1;
    Ndash = 2*N-1;
    Y(M:Mdash,N:Ndash) = X;
    X = fliplr(X);
    Y(M:Mdash,1:N-1) = X(:,1:end-1);
    X = flipud(X);
    Y(1:M-1,1:N-1) = X(1:end-1,1:end-1);
    X = fliplr(X);
    Y(1:M-1,N:Ndash) = X(1:end-1,:);
    
    if nargin > 1
        Y(~isfinite(Y)) = replaceinf;
    end
end