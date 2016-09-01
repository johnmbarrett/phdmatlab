function [times,indices] = simulatePoissonProcess(lambda,T,deltaT)
% based on https://samuelcheng.wordpress.com/2014/02/06/simulating-poisson-process-in-matlab/
    N = T/deltaT;
    
    sl = size(lambda);
    nl = numel(lambda);
    
    lambda = lambda(:)';
    
    R = rand(N,nl);
    
    events = R < repmat(lambda,N,1)*deltaT;
    
    indices = cell(nl,1);
    times = cell(nl,1);
    
    for ii = 1:nl
        indices{ii} = find(events(:,ii));
        times{ii} = indices{ii}*deltaT;
    end
    
    indices = reshape(indices,sl);
    times = reshape(times,sl);
end