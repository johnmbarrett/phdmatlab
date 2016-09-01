function I = biasCorrect(x,y,varargin)
%     assert(isequal(numel(x),size(y,1)),'There must be the same number of observations for each variable');
    
    if isempty(x) || isempty(y)
        I = 0;
        return;
    end
    
    nx = numel(x);
    ny = size(y,1);
    
    options = getopt('method=''discrete'' nreps=10 maxsplit=6 triallength=1.0 label=''ch1'' samplingfrequency=25000 stratstrat=2 singstrat=0 minembed=1 maxembed=2',varargin{:});
    
    I1 = estimateMI(x,y,options.method,varargin{:});
    
    if ~isfinite(I1)
        % splitting won't give us a more accurate estimate
        I = I1;
        return;
    end
    
    maxSplit = min(min(numel(x),numel(y)),options.maxsplit);
    n = ((0:maxSplit)/2).*((0:maxSplit)+1);
    s = zeros(n(end),1);
    
    for ii = 1:maxSplit
        s(n(ii)+1:n(ii+1)) = repmat(ii,ii,1);
    end
    
    I0 = zeros(options.nreps,1);
    Is = zeros(size(s));
    
    Is(1) = I1;
    
%     poly3 = @(a,x) a(1)*x.^2 + a(2)*x + a(3);

    badSplits = false(size(Is));
    
    for ii = 1:options.nreps
        xpermutation = randperm(nx);
        ypermutation = randperm(ny);
        
        for jj = 2:maxSplit
            k = ceil(nx/jj);
            l = ceil(ny/jj);
            
            for kk = 1:jj
                xfirst = (kk-1)*k+1;
                xlast = min(nx,kk*k);
                xindices = xpermutation(xfirst:xlast);
                
                if nx == ny % TODO : same size samples but distinct trials?
                    yindices = xindices;
                else
                    yfirst = (kk-1)*l+1;
                    ylast = min(ny,kk*l);
                    yindices = ypermutation(yfirst:ylast);
                end
                
                if isempty(xindices) || isempty(yindices)
                    badSplits(n(jj)+kk) = true;
                    continue;
                end
                
                Is(n(jj)+kk) = estimateMI(x(xindices),y(yindices,:),options.method,varargin{:});
            end
        end 
        
        % if dodgy values are introduced by the splitting, it's probably
        % just bad luck, so exclude
        p1 = polyfit(s(~badSplits & isfinite(Is)),Is(~badSplits & isfinite(Is)),2);
    %         p2 = lsqcurvefit(poly3,p1,s,Is,[-Inf -Inf 0],[Inf Inf Inf]);
    %         
    % %         if isContinuous
    %             clf;
    %             scatter(s,Is);
    %             hold on;
    %             fplot(@(x) polyval(p1,x),[0 maxSplit],'Color','r');
    %             fplot(@(x) poly3(p2,x),[0 maxSplit],'Color','g');
    % %         end

        I0(ii) = p1(3);
    end
    
    I = mean(I0);
    
%     sumN = @(x) (x/2).*(x+1);
%     n = nReps*(sumN(maxSplit)-1)+1;
%     
%     s = zeros(n,1);
%     s(1) = 1;
%     
%     I = zeros(n,1);
%     I(1) = discreteMutualInformation(x,y);
%     
%     nx = numel(x);
%     
%     for ii = 2:maxSplit
%         m = nReps*(sumN(ii-1)-1);
%         s(1+m+(1:nReps*ii)) = kron(repmat(ii,ii,1),ones(nReps,1));
%         
%         for jj = 1:nReps
%             permutation = randperm(nx);
%             k = ceil(nx/ii);
%             
%             for kk = 1:ii
%                 first = (kk-1)*k+1;
%                 last = min(nx,kk*k);
%                 indices = permutation(first:last);
%                 I(1+m+(jj-1)*ii+kk) = discreteMutualInformation(x(indices),y(indices));
%             end
%         end
%     end
%     
%     p = polyfit(s,I,2);
%     
%     I = p(3);
end

function I = estimateMI(x,y,method,varargin)
    if strncmpi(method,'f',1)
        options = getopt('mifunction=NaN',varargin{:});
        
        if ~isa(options.mifunction,'function_handle')
            error('No mutual information function supplied');
        end
        
        I = options.mifunction(x,y,varargin{:});
    elseif strncmpi(method,'m',1)
        I = binlessInfo(x,y,varargin{:});
    elseif strncmpi(method,'c',1)
        z = ~isfinite(y);
        y = num2cell(y);
        y(z) = cell(sum(z),1);
        
        I = binlessInfo(x,y,varargin{:});
    elseif strncmpi(method,'d',1)
        valid = ~isnan(x) & ~any(isnan(y),2);
        [~,~,x] = unique(x(valid));
        [~,~,y] = unique(y(valid,:),'rows'); % TODO : check the maths on this for size(y,2) > 1
        
        I = discreteMutualInformation(x,y,varargin{:});
    else
        error('Unknown method ''%s''\n',method);
    end
end