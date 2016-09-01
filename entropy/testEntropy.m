function [meanh,semh] = testEntropy(fn,maxN,repeats,h)
    hbar = zeros(maxN,repeats);
    
    for ii = 1:repeats
        for jj = 1:maxN
            x = normrnd(0,1,jj,1);
            hbar(jj,ii) = fn(x);
        end
    end
    
    if nargin > 3
        hbar = hbar-h;
    end
    
    meanh = mean(hbar,2);
    semh = std(h,[],2)/sqrt(repeats);
end

