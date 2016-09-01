function [h,mleN,ciN,mleT,ciT] = poisson_fixed_time_test(lambda,n,T)
    if nargin < 3
        T = 1;
    end
    
    if nargin < 2
        n = 1000;
    end
    
    if nargin < 1
        lambda = 1;
    end
    
    ts = poisson_fixed_time(lambda*ones(n,1),T);
    
    nts = cellfun(@numel,ts);
    
    [mleN,ciN] = mle(nts/T,'Distribution','poiss');
    
    dts = cell(size(ts));
    dts(nts > 1) = cellfun(@diff,ts(nts > 1),'UniformOutput',false);
    
    [mleT,ciT] = mle(vertcat(dts{:}),'Distribution','exp');
    mleT = 1/mleT;
    ciT = flipud(1./ciT);
    
    h = ciN(1) <= lambda && lambda <= ciN(2) && ciT(1) <= lambda && lambda <= ciT(2);
end