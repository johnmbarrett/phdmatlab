function likelihood = getKDELikelihood(response,model)
    k = size(model.prks,2);
    m = size(model.prks,3);
    
	likelihood = zeros(k,m);
    
    for ii = 1:k
        likelihood(ii,:) = interp1(model.x,model.prks(:,ii,:),response(ii));
    end
    
    % outside the bounds of the KDE
    likelihood(isnan(likelihood)) = 0;
    
    likelihood = log(likelihood);
end