function params = initDiscreteDecoder(s,R,varargin)
    m = max(s);
    n = size(R,1);
    k = size(R,2);
    binR = cell(1,k);
    vk = cell(k,1);
    
    nk = zeros(k,1);

    for ii = 1:k
        [v,~,r] = unique(R(:,ii));
        vk{ii} = v;
        nk(ii) = numel(v);
        binR{ii} = zeros(n,nk(ii),m);
        
        for jj = 1:n
            binR{ii}(jj,r(jj),s(jj)) = 1;
        end
        
        R(:,ii) = r;
    end
    
    params = struct('m',m,'nk',nk,'R',R);
    params.vk = vk;
    params.binR = binR;
end