function h = histogramEntropy(x)
    n = numel(x);
    alpha = 3;
    p = 1;
    a = (alpha-p)/(alpha*(p+2)-p^2);
    w = n^-a;
    
    objFn = @(arg) -estimateEntropy(x,arg);
    options = optimset(optimset(@fminunc),'LargeScale','off');
    
    w = fminunc(objFn,w,options);
    
    [~,h] = estimateEntropy(x,w);
end

function [Icaron,Ihat] = estimateEntropy(x,w)
	n = numel(x);
    bins = ceil(n/w);
    [f,bins] = hist(x,bins);
    Ihat = sum(plogp(f))/n - log2(n*mean(diff(bins)));
    Q = sum(f ~= 0)/n;
    Icaron = Ihat - Q;
end