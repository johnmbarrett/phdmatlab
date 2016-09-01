function I = cond2DInfoInd(s,R,varargin)
    Ps = diff(ecdf(s));
    [v,~,s] = unique(s);
    nv = numel(v);
    
    d = size(R,2);
    X = cell(d,1);
    
    for ii = 1:d
        [X{ii},~,R(:,ii)] = unique(R(:,ii));
    end
    
    nx = cellfun(@numel,X,'UniformOutput',false);
    Pr_s = cellfun(@(n) zeros(n,nv),nx,'UniformOutput',false);
    
    for ii = 1:nv
        R_s = R(s == ii,:);
        
        for jj = 1:d
            Pr_s{jj}(unique(R_s),ii) = diff(ecdf(R_s(:,jj)));
        end
    end
    
    N = prod([nx{:}]);
    PR_s = zeros(nv,nx{:});
    subs = cell(1,d);

    for hh = 1:nv
        for ii = 1:N
            [subs{:}] = ind2sub([nx{:}],ii);
            
            PR_s(hh,subs{:}) = prod(cellfun(@(p,jj) p(jj,hh),Pr_s,subs'));
        end
    end
    
    PR = permute(sum(repmat(Ps,[1 nx{:}]).*PR_s,1),[2:d+1 1]);

    I = sum(Ps.*sum(plogp(reshape(PR_s,nv,N)),2),1)-sum(plogp(PR(:)));

    options = getopt('base=2',varargin{:});

    I = I/log(options.base);
end