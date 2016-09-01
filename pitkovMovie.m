function [f,c] = pitkovMovie(N, T, dT, nu, contrast, imageSize, retinaSize)
    if nargin < 1
        N = 16;
    end
    
    if nargin < 2
        T = 256;
    end
    
    if nargin < 3
        dT = 0.015;
    end
    
    if nargin < 4
        nu = 10; % degrees visual angle per second
    end
    
    if nargin < 5
        contrast = 0.35;
    end
    
    % units don't matter so long as they're consisten
    if nargin < 6
        imageSize = 1.84;
    end
    
    if nargin < 7
        retinaSize = 17;
    end
    
    Ndash = 2*N-1;
    
    Sz = zeros(Ndash,Ndash,T);
    E = zeros(size(Sz));
    for tt = 1:T
        for ii = 0:N-1
            for jj = 0:N-1
                % ii and jj represent cycles per image
                k2 = ii^2 + jj^2;
%                 Sz(ii,jj,tt) = normrnd(0,1/k2);
%                 Sz(mod(Ndash-ii+1,Ndash)+1,mod(Ndash-jj+1,Ndash)+1,tt) = Sz(ii,jj,tt);
%                 E(ii,jj,tt) = exp(-tt/sqrt(k2));
%                 E(mod(Ndash-ii+1,Ndash)+1,mod(Ndash-jj+1,Ndash)+1,tt) = E(ii,jj,tt);
                p = normrnd(0,1/k2);
                Sz(ii+N,jj+N,tt) = p;
                Sz(N-ii,jj+N,tt) = p;
                Sz(ii+N,N-jj,tt) = p;
                Sz(N-ii,N-jj,tt) = p;
                W = atan2(imageSize,retinaSize); % visual angle subtended by image in radians
                V = 180*W/pi; % visual angle in degrees
                % therefore sqrt(k2)/V gives cycles per degree
                e = exp(-dT*tt/(nu*sqrt(k2)/V));
                E(ii+N,jj+N,tt) = e;
                E(N-ii,jj+N,tt) = e;
                E(ii+N,N-jj,tt) = e;
                E(N-ii,N-jj,tt) = e;
            end
        end
        Sz(N,N,tt) = 0;
        Sz(:,:,tt) = ifftshift(squeeze(Sz(:,:,tt)));
        E(:,:,tt) = ifftshift(squeeze(E(:,:,tt)));
        
%         Stt = S(1:N,1:N,tt);
%         St = reshape(Stt,numel(Stt),1);
%         
%         mu = mean(St);
%         sigma = std(St);
%         
%         if mu == 0
%             % panic
%             tt = tt-1; %#ok<FXSET>
%             continue;
%         end
%         
%         A = sigma/mu*contrast;
%         
%         S(:,:,tt) = A*S(:,:,tt);
    end
    
    S = zeros(size(Sz));
    
    for ii = 1:Ndash
        for jj = 1:Ndash
            S(ii,jj,:) = conv(squeeze(E(ii,jj,:)),squeeze(Sz(ii,jj,:)),'same');
        end
    end
    
    1+1;
    
    s = zeros(size(S));
    for ii = 1:T
        s(:,:,ii) = ifft2(S(:,:,ii));
        
        1+1;
    end
    
    c = zeros(N,N,T);
    h = figure;
    for ii = 1:T
        d = s(1:N,1:N,ii);
        c(:,:,ii) = d;
        dmin = min(min(d));
        dmax = max(max(d));
        image(255*(d-dmin)/(dmax-dmin));
        colormap('gray');
        f(ii) = getframe;
    end
    
    close(h);
end