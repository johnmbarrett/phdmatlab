function I = continuousMutualInformation(X,k,estimator)
    if nargin < 3
        estimator = 1;
    elseif ~isnumeric(estimator) || ~isscalar(estimator) || estimator < 1 || estimator > 2
        error('Estimator must be 1 or 2');
    end
    
    if nargin < 2
        k = 1;
    end
    
    if isscalar(X) || isvector(X)
        error('Mutual information between a random variable and itself is just its entropy, which this estimator isn''t designed for.');
    elseif ismatrix(X)
        d = 1;
        [N,m] = size(X);
    elseif ndims(X) > 3
        error('what ar ey ou trying to do just no????');
    else
        [N,d,m] = size(X);
    end
    
    assert(k < N,'kNN parameter must be stricly less than the number of samples.');
    
    X = reshape(X,[N 1 d m]);
    
    rho = sort(squeeze(sum((repmat(X,1,N,1)-repmat(permute(X,[2 1 3 4]),N,1,1)).^2,3)),2);
    
    switch estimator
        case 1
            epsilon = squeeze(max(squeeze(rho(:,k+1,:)),[],2));
            n = squeeze(sum(rho(:,2:end,:) < repmat(epsilon,1,N-1,m),2));

            I = psi(k) + (m-1)*psi(N) - mean(sum(psi(n+1),2),1);
        case 2
            epsilon = rho(:,k+1,:);
            n = squeeze(sum(rho(:,2:end,:) <= repmat(epsilon,[1 N-1 1]),2));
            
            I = psi(k) - (m-1)*k + (m-1)*psi(N) - mean(sum(psi(n),2),1);
    end
end