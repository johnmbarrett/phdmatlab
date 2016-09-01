function I = contrastGratingsMutualInformationDoubleNeuron(sigma,mu,lambda)
%     lmax = max([mu; lambda(:)]);
%     xmax = poissinv(1-eps,lmax);
%     pmin = poisspdf(xmax,lmax);
%     
%     tic;
%     while pmin > 0
%         xmax = xmax*2;
%         pmin = poisspdf(xmax,lmax);
%     end
%     toc;
%     
%     tic;
    x = (0:14085)'; % poisspdf(x(end),L) guaranteed to be 0 for L <= 10000
    
    mpdf = poisspdf(x,mu+sigma);
    gpdf = arrayfun(@(l) poisspdf(x,l+sigma),lambda,'UniformOutput',false);
    
    Hs = cellfun(@(pdf) -mean(sum([plogp(pdf) plogp(mpdf)]),2),gpdf);
    H = cellfun(@(pdf) -sum(plogp(mean([pdf mpdf],2))),gpdf);

    I = H - Hs;
%     toc;
end