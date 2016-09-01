function dI = deltaI(s,R,varargin)
    [v,~,s] = unique(s);
    ns = accumarray(s,1);
    Ps = ns/numel(s);
    
    d = size(R,2);
    X = cell(d,1);
    
    for ii = 1:d
        [X{ii},~,R(:,ii)] = unique(R(:,ii));
    end
    
    nv = numel(v);
    nx = cellfun(@numel,X)'; %,'UniformOutput',false);
    
    PR_s = zeros([nv nx]);
    Pr_s = arrayfun(@(n) zeros(nv,n),nx,'UniformOutput',false);
    colons = repmat({':'},1,d);
    
    for ii = 1:nv
        R_s = R(s == ii,:);
        ns = size(R_s,1);
        PR_s(ii,colons{:}) = accumarray(R_s,1,nx)/ns;
        
        for jj = 1:d
            Pr_s{jj}(ii,:) = accumarray(R_s(:,jj),1,[nx(jj) 1])/ns;
        end
    end
    
%     tic;
    QR_s = ones(size(PR_s));
    
    for ii = 1:d
        reps = [1 nx];
        reps(ii+1) = 1;
        shape = [nv ones(1,d)];
        shape(ii+1) = nx(ii);
        QR_s = QR_s.*repmat(reshape(Pr_s{ii},shape),reps);
    end
%     toc;
    
%     tic;
%     M = prod(nx);
%     QR_s = ones(nv,M);
%     cx = num2cell(nx);
%     subs = cell(1,d);
%     
%     for ii = 1:M
%         [subs{:}] = ind2sub([cx{:}],ii);
%         
%         for jj = 1:d
%             sub = subs{jj};
%             QR_s(:,ii) = QR_s(:,ii).*Pr_s{jj}(:,sub);
%         end
%     end
%     toc;
    
    M = prod(nx);
    PR_s = reshape(PR_s,nv,M);
    QR_s = reshape(QR_s,nv,M);
    
    PR = sum(PR_s.*repmat(Ps,1,M))';
    
    % won't strictly be true due to floating point errors and in any case
    % omitted for speed, but a useful sanity check
%     assert(isequal(PR,reshape(accumarray(R,1,nx)/size(R,1),size(PR))));
    
    QR = sum(QR_s.*repmat(Ps,1,M))';
    
    dI = sum(Ps.*sum(plogq(PR_s,PR_s./QR_s),2),1)-sum(plogq(PR,PR./QR),1);
    
    options = getopt('base=2',varargin{:});

    dI = dI/log(options.base);
    
    return;
    
    Ps = diff(ecdf(s));
    [v,~,s] = unique(s);
    nv = numel(v);
    
    d = size(R,2);
    X = cell(d,1);
    
    for ii = 1:d
        [X{ii},~,R(:,ii)] = unique(R(:,ii));
    end
    
    nx = cellfun(@numel,X,'UniformOutput',false);
    
    N = prod([nx{:}]);
    disp(N);
    PR_s = zeros(nv,N);
    QR_s = zeros(nv,nx{:});
    subs = cell(1,d);
    
    for ii = 1:nv
%         tic;
        R_s = R(s == ii,:);
        ns = size(R_s,1);
%         pr_s = cellfun(@(n) zeros(n,1),nx,'UniformOutput',false);
%         
%         for jj = 1:d
%             pr_s{jj}(unique(R_s(:,jj))) = diff(ecdf(R_s(:,jj)));
%         end
        QR_sfun = @(kk,val) sum(R_s(:,kk) == val)/ns;
%         toc;
        
%         tic;
        for jj = 1:N
            tic;
            [subs{:}] = ind2sub([nx{:}],jj);
            toc;
            tic;
            PR_s(ii,jj) = sum(ismember(R_s,[subs{:}],'rows'))/ns;
            toc;
            tic;
            QR_s(ii,jj) = prod(cellfun(QR_sfun,num2cell(1:d),subs));
            toc;
%             QR_s(ii,subs{:}) = prod(cellfun(@(p,kk) p(kk),pr_s,subs'));
        end
%         toc;
    end
    
    QR_s = reshape(QR_s,nv,N);
    QR = sum(QR_s.*repmat(Ps,1,N))';
    PR = sum(PR_s.*repmat(Ps,1,N))';
    
    dI = sum(Ps.*sum(plogq(PR_s,PR_s./QR_s),2),1)-sum(plogq(PR,PR./QR),1);
    
    options = getopt('base=2',varargin{:});

    dI = dI/log(options.base);
end