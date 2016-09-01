function likelihood = getJacobsLikelihood(response,model)
    likelihood = -model.Vs; % log(exp(Vs)) == Vs
    
    bw = model.bw;
    response = cellfun(@(r) ceil(r/bw),response,'UniformOutput',false);
    
    vs = model.vs;
    nr = size(vs,2);
    
    for ii = 1:nr
        likelihood(ii,:) = likelihood(ii,:) + permute(sum(log(vs(response{ii},ii,:)*bw),1),[2 3 1]);
    end
end