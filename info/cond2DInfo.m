function I = cond2DInfo(s,R,varargin)
    Ps = diff(ecdf(s));

    [~,~,r] = unique(R,'rows');
    PR = diff(ecdf(r));

    [v,~,s] = unique(s);
    nv = numel(v);

    PR_s = zeros(size(PR,1),nv);

    for hh = 1:nv
        r_s = r(s == hh);
        PR_s(unique(r_s),hh) = diff(ecdf(r_s));
    end

    I = sum(Ps'.*sum(plogp(PR_s),1),2)-sum(plogp(PR));

    options = getopt('base=2',varargin{:});

    I = I/log(options.base);
end