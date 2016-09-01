function likelihood = getMBGLikelihood(response,model)
    k = size(model.params,1);
    m = size(model.params,2);
    
    likelihood = zeros(k,m);
    
    for ii = 1:k
        if isinf(response(ii))
            likelihood(ii,:) = model.params(ii,:,1);
            continue;
        end
        
        likelihood(ii,:) = (1-model.params(ii,:,1)).*normpdf(response(ii),model.params(ii,:,2),model.params(ii,:,3));
        
        isDegenerate = model.params(ii,:,3) == 0;
        
        % MATLAB gives normpdf(x,m,0) == NaN for all x,m - but we can do
        % better than that. A 0 std normal distribution is a Dirac delta,
        % so p(x == m) = 1 and p(x == n) = 0 for all n ~= m.
        likelihood(ii,isDegenerate) = response(ii) == model.params(ii,isDegenerate,2);
    end
    
    % NaNs will happen when the test response is finite but there were only
    % Infs in the training set (because (1-1)*normpdf(x,NaN,NaN) = NaN).
    likelihood(isnan(likelihood)) = 0;
    
    likelihood = log(likelihood);
end