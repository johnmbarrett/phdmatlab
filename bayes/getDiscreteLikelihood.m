function psrk = getDiscreteLikelihood(response,pdfs)
    prks = pdfs.prks;    

    m = size(prks{1},2);
    k = numel(response);
    psrk = zeros(k,m);
        
    for ii = 1:k
        idx = response(ii) == pdfs.vk{ii};
        
        if ~any(idx)
            psrk(ii,:) = 0;
            continue;
        end
        
        psrk(ii,:) = prks{ii}(idx,:);
    end
    
    psrk = log(psrk);
end