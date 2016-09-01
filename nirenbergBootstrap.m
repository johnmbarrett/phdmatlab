function [k,pks] = nirenbergBootstrap(data,interval,bw,tol,N,alpha)
    if numel(data) == 1
        k = 1;
        return;
    end
    
    if nargin < 7
        alpha = 0.05;
    end
    
    if nargin < 6
        N = 500;
    end

    if nargin < 5
        tol = 0.001;
    end
    
    if nargin < 3
        bw = 0.001;
    end
    
    if nargin < 2
        omega = @(z,h) 1;
    else
        a = min(interval);
        b = max(interval);
        
        if a == b
            error('diff(range) must not be zero');
        end
        
        omega = @(z,h) normcdf(b,z,h) - normcdf(a,z,h);
    end
    
    k = 0;
    pk = 0;
    pks = [];
    
    while pk < alpha
        k = k + 1;
        
        h1 = findSmallestH(data,k,omega,bw,tol);

        [p,X] = empPDF(data,h1,omega,bw);
        P = cumsum(p);

        hs = zeros(N,1);
        for ii = 1:N
            x = rand(size(data));
            q = interp1(P,X,x);
            hs(ii) = findSmallestH(q,k,omega,bw,tol);
        end

        pk = sum(hs > h1)/N;
        pks = [pks; pk]; %#ok<AGROW>
    end
end

function hmax = findSmallestH(data,k,omega,bw,tol)
    hnew = 1;
    knew = countModes(data,hnew,omega,bw);
    
    if knew <= k
        while knew <= k
            hold = hnew;
            hnew = hnew/2;
            knew = countModes(data,hnew,omega,bw);
        end
    else
        while knew > k
            hold = hnew;
            hnew = hnew*2;
            knew = countModes(data,hnew,omega,bw);
        end
    end
    
    hmin = min([hold hnew]);
    hmax = max([hold hnew]);
    
    while abs(hmax-hmin) > tol
        hmid = mean([hmin hmax]);
        
        kmid = countModes(data,hmid,omega,bw);
        
        if kmid <= k
            hmax = hmid;
        else
            hmin = hmid;
        end
    end
                
%         kmax = countModes(data,hmin);
%         
%         if kmax == 1
%             hmin = hmin/2;
%             continue;
%         end
%         
%         kmin = countModes(data,hnew,omega,bw);
%         
%         if kmin == 1
%             hnew = mean([hmin hnew]);
%         else
%             hnew = 2*hnew-hmin;
%         end
end

function k = countModes(data,h,omega,bw)
    if h == 0
        k = numel(data);
        return;
    end
    
    [P,X,F] = empPDF(data,h,omega,bw);
    
    plotyy(X((1:numel(F))+ceil((numel(P)-numel(F))/2)),F,X,P);
    
    k = numel(findpeaks(P));
end

function [P,X,F] = empPDF(data,h,omega,bw)
    X = (floor(min(data)/bw)*bw:bw:ceil(max(data)/bw)*bw)';
    F = zeros(size(X));

    for ii = 1:numel(data)
        F(floor((data(ii)-min(data))/bw)+1) = omega(data(ii),h);
    end
    
    kernel = normpdf((-ceil(5*h/bw):ceil(5*h/bw))'*bw,0,h);
    
    assert(mod(numel(kernel),2) == 1,'Kernel must have an odd number of elements');
    n = floor(numel(kernel)/2);
    
    P = conv(F,kernel);
    P = P/sum(P);
    
    X = [X(1)-(n:-1:1)'*bw; X; X(end)+(1:n)'*bw];
    
    assert(numel(X) == numel(P),'PDF and domain must have same number of elements');
end