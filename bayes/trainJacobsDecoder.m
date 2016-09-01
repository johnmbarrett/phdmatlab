function modelFun = getJacobsDecoderTrainingFun(bw,sd)
    modelFun = @(idx,s,R) trainJacobsDecoder(idx,s,R,bw,sd);
end

function model = trainJacobsDecoder(trainIndices,s,R,bw,sd)
    ns = max(s);
    nr = size(R,2);
    
    T = ceil(1000*max(max(cellfun(@(t) tertiaryop(isempty(t),0,max(t)),R(trainIndices,:)))));
    hists = zeros(T,ns,nr);

    for ii = 1:ns
        stimulusTrials = find(trainIndices & s == ii);
        nt = numel(stimulusTrials);
        
        for jj = 1:nr
            for kk = 1:nt
                row = ceil(1000*R{stimulusTrials(kk),jj});
                hists(row,ii,jj) = hists(row,ii,jj) + 1;
            end
        end
    end

    %%

    ifrs = cellfun(@(A) zeros(size(A)),hists,'UniformOutput',false);
    kernel = normpdf(-125:125,0,25);
    kernel = kernel/sum(kernel);

    %%

    for ii = 1:2
        for jj = 1:nCells
            for kk = 1:4
                ifrs{ii}(:,kk,jj) = conv
end