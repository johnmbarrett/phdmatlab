function likelihood = getMBKLikelihood(response,model)
    k = size(model.qrks,1);
    
    likelihood = getKDELikelihood(response,model);
    
    for ii = 1:k
        if isinf(response(ii))
            likelihood(ii,:) = log(model.qrks(ii,:));
            continue;
        end
        
        likelihood(ii,:) = log(1-model.qrks(ii,:))+likelihood(ii,:);
    end
end