function sample = bootstrap(data,nTrials,iter,useMedian)
    avg = @mean;
    
    if nargin > 3 && useMedian
        avg = @median;
    end
    
    if nargin < 3
        iter = 10000;
    end
    
    if nargin < 2
        nTrials = 10;
    end
    
    if nargin < 1
        error('You must provide a sample to bootstrap from');
    end
    
    n = numel(data);
    data = data(:);
    
    nConditions = n/nTrials;
    
    if mod(nConditions,1) ~= 0
        error('Each condition must have the same number of trials');
    end
    
    sample = zeros(iter,nConditions);
    
    for qq = 1:iter
        sample(qq,:) = avg(reshape(data(randperm(n)),nTrials,nConditions));
    end
    
    m = nConditions*iter;
    sample = reshape(sample,m,1);
end