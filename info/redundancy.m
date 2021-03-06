function R = redundancy(s,r)
    assert(size(s,1) == size(r,1),'There must be the same number of stimulus and response observations.');
    
    N = size(s,1);
    
    nSources = size(r,2);
    
    stimuli = unique(s);
    nStimuli = numel(unique(s));
    
    stimulusPDF = zeros(nStimuli,1);
    specificInformation = zeros(nStimuli,nSources);
    
    for ii = 1:nStimuli
        si = s == stimuli(ii);
        stimulusPDF(ii) = sum(si)/N;
        
        for jj = 1:nSources
            z = r(:,jj);
            nz = size(z,1);
            
            zcond = r(si,jj);
            nzcond = size(zcond,1);
            
            pInf = sum(isinf(z))/nz;
            pInfCond = sum(isinf(zcond))/nzcond;
            
            Dx = pInf*log(pInf/pInfCond) + (1-pInf)*log((1-pInf)/(1-pInfCond));
            Dy = continuousKLDivergence(z(isfinite(z)),zcond(isfinite(zcond)),1,'both');
            
            specificInformation(ii,jj) = Dx + (1-pInf)*Dy;
        end
    end
        
    R = sum(stimulusPDF.*min(specificInformation,[],2));
end