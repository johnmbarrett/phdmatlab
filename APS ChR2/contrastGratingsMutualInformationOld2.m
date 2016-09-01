function I = contrastGratingsMutualInformation(sigma,p,c)
    if nargin < 3
        c = 0.1:0.1:0.6;
    end
    
    lambdaFun = @(x) p(1).*normcdf(x,p(2),p(3)); %@(x) p(1)./(1+exp(-(x-p(2))/p(3)));
    mu = lambdaFun(1); %polyval(p,m);
    lambda = lambdaFun([1+c;1-c]); %max(0,polyval(p,[m*(1+c);m*(1-c)]));

    %%

    I = zeros(2,numel(c));
    n = 0:170;

    for ii = 1:numel(c)
        for jj = 1:2
             mpdf = @(x) mean([poisspdf(x,sigma+mu); poisspdf(x,sigma+lambda(jj,ii))]);

             I(jj,ii) = -sum(plogp(mpdf(n)))                ...
                 + sum(plogp(poisspdf(n,sigma+lambda(jj,ii))))/2  ...
                 + sum(plogp(poisspdf(n,sigma+mu)))/2;
        end
    end

    I = sum(I)';
end