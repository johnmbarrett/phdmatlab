function I = contrastGratingsMutualInformation(sigma,mu,p)
%     sigma = 1;
    c = polyval(p,(1:6));
%     mu = 3;
    lambda = [mu+c;max(0,mu-c)];

    %%

    I = zeros(2,6);
    n = 0:170;

    for ii = 1:6
        for jj = 1:2
             mpdf = @(x) mean([poisspdf(x,sigma+mu); poisspdf(x,sigma+lambda(jj,ii))]);

             I(jj,ii) = -sum(plogp(mpdf(n)))                ...
                 + sum(plogp(poisspdf(n,sigma+lambda(jj,ii))))/2  ...
                 + sum(plogp(poisspdf(n,sigma+mu)))/2;
        end
    end

    I = sum(I);
end