function [f,c] = pitkowMovie(N, T, dT, nu, contrast, imageSize, retinaSize)
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
    
    W = atan2(imageSize,retinaSize); % visual angle subtended by image in radians
    V = 180*W/pi; % visual angle in degrees
    [I,J] = meshgrid((0:N-1).^2);
    sigmas = 1./(I+J);
                
    for tt = 1:T
        Sz(:,:,tt) = ifftshift(fft2reflect(normrnd(0,sigmas),0));
        E(:,:,tt) = ifftshift(fft2reflect(exp(-dT*tt./(nu*sqrt(I+J)/V))));
        surf(squeeze(E(:,:,1)));
        
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
        colormap(gray(256));
        f(ii) = getframe;
    end
    
    close(h);
end