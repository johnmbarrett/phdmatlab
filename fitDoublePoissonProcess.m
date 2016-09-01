function [isiThresh,lambdaShort,lambdaLong,mse] = fitDoublePoissonProcess(isis)
    function mse = modelFn(isiThresh,isis)
        shortISIs = isis(isis < isiThresh);
        longISIs = isis(isis >= isiThresh)-max(shortISIs);
        
        mus = expfit(shortISIs);
        mul = expfit(longISIs);
        
        [fs,xs] = ecdf(shortISIs);
        [fl,xl] = ecdf(longISIs);
        
        mses = mean((fs-expcdf(xs,mus)).^2);
        msel = mean((fl-expcdf(xl,mul)).^2);
        
        mse = mses+msel;
    end

    sortedISIs = sort(isis);
    isiThresh = isis(ceil(numel(sortedISIs)/2));
    
    [isiThresh,mse] = fmincon(@(x) modelFn(x,isis),isiThresh,[],[],[],[],sortedISIs(2),max(isis));
    
    lambdaShort = 1/expfit(isis(isis < isiThresh));
    lambdaLong = 1/expfit(isis(isis >= isiThresh));
end