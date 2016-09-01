function I = contrastGratingsMutualInformation(sigma,p,c,mu,lambda,xmax)
    if nargin < 6
        xmax = 14085;
    end
    
    if nargin >= 5
        nn = size(lambda,2);
        assert(nn == 2,'Only two neurons are supported at present'); % TODO : n > 2
        nc = size(lambda,1);
        
        if numel(sigma) == 1
            sigma = repmat(sigma,1,nn);
        else
            assert(numel(sigma) == nn,'Sigma must either be a scalar or a vector with one value per neuron');
        end
        
        if numel(mu) == 1
            mu = repmat(mu,1,nn);
        else
            assert(numel(mu) == nn,'Sigma must either be a scalar or a vector with one value per neuron');
        end
    else
        if nargin < 3
            c = 0.1:0.1:0.6;
        end

        lambdaFun = @(x) p(1).*normcdf(x,p(2),p(3)); %@(x) p(1)./(1+exp(-(x-p(2))/p(3)));
        mu = lambdaFun(1); %polyval(p,m);
        lambda = lambdaFun([1+c;1-c])'; %max(0,polyval(p,[m*(1+c);m*(1-c)]));
        nc = numel(c);
    end

    %%
    
    Hs = zeros(nc,1);
    H = zeros(nc,1);
    n = 0:xmax;
    
    tic;
    mpdf1 = poisspdf(n',mu(1)+sigma(1));
    
    if mu(1)+sigma(1) == mu(2)+sigma(2)
        mpdf2 = mpdf1';
    else
        mpdf2 = poisspdf(n,mu(2)+sigma(2));
    end
    
    mjpdf = repmat(mpdf1,1,xmax+1).*repmat(mpdf2,xmax+1,1);
    
    for ii = 1:nc
        gpdf1 = poisspdf(n',lambda(ii,1)+sigma(1));
        
        if lambda(ii,1)+sigma(1) == lambda(ii,2)+sigma(2)
            gpdf2 = gpdf1';
        else
            gpdf2 = poisspdf(n,lambda(ii,2)+sigma(2));
        end
        
        gjpdf = repmat(gpdf1,1,xmax+1).*repmat(gpdf2,xmax+1,1);
        
        Hs(ii) = -mean(sum([plogp(gjpdf(:)) plogp(mjpdf(:))]),2);
        H(ii) = -sum(plogp(mean([gjpdf(:) mjpdf(:)],2)));
    end
    toc;

    I = H - Hs;
end