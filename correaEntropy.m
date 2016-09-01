function Cmn = correaEntropy(x,m,unbiased)
    n = numel(x);
    
    if unbiased
        CmnHat = correaEntropy(x,m,false);
        
        CmnHatI = zeros(n,1);
        
        for ii = 1:n
            CmnHatI(ii) = correaEntropy(x(setdiff(1:n,ii)),m,false);
        end
        
        CmnHatSigma = mean(CmnHatI);
        bias = (n-1)*(CmnHatSigma-CmnHat);
        Cmn = CmnHat - bias;
        
        return;
    end
    
    y = sort(x);
    y = [y(1)*ones(m,1); y; y(end)*ones(m,1)];
    b = zeros(size(x));
    
    for ii = 1:n
        yj = y((0:2*m)+ii);
        yibar = mean(yj);
        z = (yj-yibar);
        b(ii) = sum(z.*(-m:m)')./(n*sum(z.^2));
    end
    
    Cmn = -sum(log(b))/n;
end