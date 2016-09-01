function [s,v,u] = izhNeuron(v0,tmax,tI,I,a,b,c,d)
    % based on code in Izhikevich (2003), IEEE Trans Neural Network
    if nargin < 8
        d = 8;
    end
    
    if nargin < 7
        c = -65;
    end
    
    if nargin < 6
        b = 0.2;
    end
    
    if nargin < 5
        a = 0.02;
    end
    
    if nargin < 2
        tmax = 1000;
    end
    
    if nargin < 4
        tI = 0:ceil(tmax);
        I = zeros(ceil(tmax)+1,1);
    end
    
    if nargin < 1
        v0 = -65;
    end
    
    ts = (0:ceil(tmax))';
    nt = numel(ts);
    nI = size(I,2);
    v = zeros(nt,nI);
    u = zeros(nt,nI);
    s = zeros(nt,nI);
    
    doInterp = size(I,1) ~= nt;
    
    v(1) = v0;
    u(1) = b*v0;
    
    for tt = 1:nt-1
%         tic;
        if doInterp
            It = interp1(tI,I,ts(tt));
        else
            It = I(tt,:);
        end
%         toc;
        
%         tic;
        isSpiking = v(tt,:) >= 30;
        s(tt,isSpiking) = 1;
        v(tt,isSpiking) = c;
        u(tt,isSpiking) = u(tt,isSpiking) + d;
%         toc;
        
%         tic;
        % step size 0.5 ms for v for numerical stability
        vmid      = v(tt,:) + 0.5*(0.04*v(tt,:).^2 + 5*v(tt,:) + 140 - u(tt,:) + It);
        v(tt+1,:) = vmid    + 0.5*(0.04*vmid.^2    + 5*vmid    + 140 - u(tt,:) + It);
%         toc;
        
%         tic;
        u(tt+1,:) = u(tt,:) + a.*(b.*v(tt+1,:) - u(tt,:));
%         toc;
    end;
end

    