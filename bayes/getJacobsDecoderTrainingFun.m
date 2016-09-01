function [modelFun,initFun] = getJacobsDecoderTrainingFun(bw,sd,T,densityEstimationMethod)
    if nargin < 4
        densityEstimationMethod = 'kde';
    end
    
    if nargin < 3
        T = NaN;
    end
    
    initFun = @(varargin) struct('bw',bw,'sd',sd/bw,'T',T,'density_method',densityEstimationMethod);
    modelFun = @trainJacobsDecoder;
end

function model = trainJacobsDecoder(trainIndices,s,R,params)
    ns = max(s);
    nr = size(R,2);
    
    bw = params.bw;
    
    if isscalar(params.T) && isnumeric(params.T) && isfinite(params.T) && params.T > 0
        T = ceil(params.T/bw);
    else
        T = ceil(max(max(cellfun(@(t) tertiaryop(isempty(t),0,max(t)),R)))/bw);
    end
    
    hists = zeros(T,nr,ns);
    nts = zeros(ns,1);

    for ii = 1:ns
        stimulusTrials = trainIndices(s(trainIndices) == ii);
        nt = numel(stimulusTrials);
        nts(ii) = nt;
        
        for jj = 1:nr
            for kk = 1:nt
                row = ceil(1000*R{stimulusTrials(kk),jj});
                hists(row,jj,ii) = hists(row,jj,ii) + 1;
            end
        end
    end

    %%
    
    ifrs = zeros(size(hists));
        sd = params.sd;
        kernel = normpdf(-5*sd:sd*5,0,sd);
        kernel = kernel/sum(kernel);

    if strcmp(params.density_method,'kde')

        for ii = 1:ns
            for jj = 1:nr
                ifrs(:,jj,ii) = conv(hists(:,jj,ii)/nts(ii),kernel,'same');
            end
        end
    elseif strcmp(params.density_method,'bars')
        for ii = 1:ns
            for jj = 1:nr
                try
                    bars = barsP(hists(:,jj,ii),[bw T*bw],nts(ii));
                    ifrs(:,jj,ii) = bars.mean/nts(ii); 
                catch err
                    logMatlabError(err);
                    warning('BARS failed, failing back to KDE');
                    ifrs(:,jj,ii) = conv(hists(:,jj,ii)/nts(ii),kernel,'same');
                end
            end
        end
    else
        error('Unknown spike density function estimation method');
    end
    
    model = struct('bw',bw,'vs',ifrs,'Vs',permute(trapz((1:T)*bw,ifrs),[2 3 1]));
end