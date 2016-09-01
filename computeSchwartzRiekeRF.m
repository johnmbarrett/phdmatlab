function computeSchwartzRiekeRF(recording,varargin)
    channelInfo = getMetaInfo(recording.index,'spont',true);
    channelIndices = parseChannelOptions(recording.index,varargin{:});
    fileDir = getAnalysisOutputDir(recording);
    
    options = optimset(optimset(@fmincon),'Algorithm','sqp','Display','off');
    
    fig = figure;
    set(fig,'Position',[0 0 1600 900],'PaperPositionMode', 'auto')
%     for ii = find(channelIndices == find(strcmp('13',{channelInfo.label})));
%     for ii = (find(channelIndices == find(strcmp('67',{channelInfo.label})))):numel(channelIndices);
    for ii = 1:numel(channelIndices)
        channelIndex = channelIndices(ii);
        channelLabel = channelInfo(channelIndex).label;
        
        stas = dir([fileDir '\sta_ ' recording.dataFile '_channel_' channelLabel '_cluster_*.gif']);
        stas = {stas.name};
        
%         escape = false;
%         for jj = 4
        for jj = 1:numel(stas)
            fprintf('Fitting RF for channel %s cluster #%d\n',channelLabel,jj);
            tic;
            sta = [fileDir '\' stas{jj}];
            
            pixels = squeeze(double(imread(sta,'frames','all')));
            pixels = pixels(1:2:end,1:2:end,:); % TODO : read pixel size from stim file
            Z = (pixels - mean(reshape(pixels,numel(pixels),1)))/std(reshape(pixels,numel(pixels),1));
            
            nX = size(pixels,2);
            nY = size(pixels,1);
            nT = size(pixels,3);
                
            [Y,X] = ndgrid(1:nY,1:nX);
            T = 1:nT;
            
%             temporalRF = ones(size(pixels,3),1);
            
%             for nn = 1:2
%                 spatialProjection = pixels.*(repmat(reshape(temporalRF,[1 1 numel(temporalRF)]),[nY nX 1]));
%                 spatialMax = max(spatialProjection,[],3);
%                 spatialMin = min(spatialProjection,[],3);
%                 
%                 imax = find(spatialMax == max(max(spatialMax)));
%                 imin = find(spatialMin == min(min(spatialMin)));
                oldZ = Z;
                
                for tt = 1:nT
                    Z(:,:,tt) = conv2(squeeze(oldZ(:,:,tt)),[1 2 1; 2 4 2; 1 2 1]/16,'same');
                end
                
                imax = find(Z == max(max(max(Z))));
                imin = find(Z == min(min(min(Z))));
                iextreme = find(abs(Z) == max(max(max(abs(Z)))));
%                 [ymax,xmax] = ind2sub(size(spatialMax),imax);
%                 [ymin,xmin] = ind2sub(size(spatialMin),imin);
                [~,~,tmax] = ind2sub(size(pixels),imax);
                [~,~,tmin] = ind2sub(size(pixels),imin);
                [yextreme,xextreme,textreme] = ind2sub(size(pixels),iextreme);
                
                Z = oldZ;
                
                fprintf('%d candidate minima, %d candidate maxima and %d candidate extrema found for channel %s cluster #%d\n',numel(tmax),numel(tmin),numel(textreme),channelLabel,jj);
            
                [~,~,~,fileSuffix] = regexp(sta,'_channel_[0-9]+_cluster_[0-9]+.*');
                saveFile = [fileDir '\srrf_' recording.dataFile fileSuffix{1}(1:end-4)];
                
                % no extrema is usually a result of NaNs, which shouldn't
                % happen, the second condition is an heuristic: this never 
                % seems to happen for valid RFs
                if any([numel(tmax) numel(tmin) numel(textreme)] == 0) || numel(tmax)*numel(tmin) > 9
%                 if max([numel(tmax); numel(tmin); numel(textreme)]) >= max(nX,nY)
                    % the equivalent of an entire line is all black or all 
                    % white; which is extremely unlikely in a valid RF
                    significant = false; %#ok<NASGU>
                    badSTA = true; %#ok<NASGU>
                    save(saveFile,'significant','badSTA');
                    fprintf('Aborted fitting RF for channel %s cluster #%d after %f seconds',channelLabel,jj,toc);
                    continue;
                end
                
%                 dmax = sqrt((yextreme-ymax).^2 + (xextreme-xmax).^2);
%                 tmax = tmax(dmax == min(dmax));
%                 
%                 dmin = sqrt((yextreme-ymin).^2 + (xextreme-xmin).^2);
%                 tmin = tmin(dmin == min(dmin));

%                 maxZ = max(Z(yextreme,xextreme,:));
%                 minZ = min(Z(yextreme,xextreme,:));
                maxZ = max(max(max(Z)));
                minZ = min(min(min(Z)));
%                 dog = ...
%                     @(x,y,tmax,sdxmax,sdymax,sdtmax,phimax, ...
%                     tmin,sdxmin,sdymin,sdtmin,phimin) ...
%                     maxZ*gauss3d(X-x,Y-y,T-tmax,sdxmax,sdymax,sdtmax,phimax,true) + ...
%                     minZ*gauss3d(X-x,Y-y,T-tmin,sdxmin,sdymin,sdtmin,phimin,true);
%                 
%                 objFn = @(args) mean(mean(mean((Z-dog(args(1),args(2),args(3),args(4),args(5),args(6),args(7),args(8),args(9),args(10),args(11),args(12))).^2)));
%                 argmin = fmincon(objFn, ...
%                     [xextreme(1) yextreme(1) tmax(1) 1 1 1 0 ...
%                     tmin(1) 1 1 1 0],[],[],[],[], ...
%                     [1 1 1 1 1 1 -pi 1 1 1 1 -pi], ...
%                     [nX nY nT nX nY nT pi nT nX nY nT pi],[],options);

%                 objFnFn = @(k) @(args) mean(mean(mean((Z-k*gauss3d(X-args(1),Y-args(2),T-args(3),args(4),args(5),args(6),args(7),true)).^2)));
                
%                 xStarts = [1 nX/2 nX];
%                 yStarts = [1 nY/2 nY];
%                 tStarts = [1 nT/2 nT];
%                 
%                 onargmins = zeros(3,7,3,3);
%                 fvals = zeros(3,3,3);
%                 
%                 for kk = 1:3
%                     for ll = 1:3
%                         for mm = 1:3
%                             [onargmins(kk,:,ll,mm),fvals(kk,ll,mm)] = fmincon(objFnFn(Z,X,Y,T,maxZ,true), ...
%                                     [xStarts(kk) yStarts(ll) tStarts(mm) 1 1 1 0],[],[],[],[], ...
%                                     [1 1 1 1 1 1 -pi],[nX nY nT nX nY nT pi],[],options);
%                         end
%                     end
%                 end
                onargmins = zeros(numel(textreme),7,numel(tmax));
                mses = zeros(numel(textreme),numel(tmax));
                
                for kk = 1:numel(textreme)
                    for ll = 1:numel(tmax)
%                         [onargmins(kk,:,ll),mses(kk,ll)] = fmincon(objFnFn(maxZ), ...
                        [onargmins(kk,:,ll),mses(kk,ll)] = fmincon(objFnFn(Z,X,Y,T,maxZ,true), ...
                            [xextreme(kk) yextreme(kk) tmax(ll) 1 1 1 0],[],[],[],[], ...
                            [1 1 1 1 1 1 -pi],[nX nY nT nX nY nT pi],[],options);
                    end
                end
                
                iminmse = find(mses == min(min(mses)));
                [kk,ll] = ind2sub(size(mses),iminmse);
                onargmin = squeeze(onargmins(kk,:,ll));
                
                offargmins = zeros(numel(textreme),7,numel(tmin));
                mses = zeros(numel(textreme),numel(tmin));
                
                for kk = 1:numel(textreme)
                    for ll = 1:numel(tmin)
%                         [offargmins(kk,:,ll),mses(kk,ll)] = fmincon(objFnFn(minZ), ...
                        [offargmins(kk,:,ll),mses(kk,ll)] = fmincon(objFnFn(Z,X,Y,T,minZ,true), ...
                            [xextreme(kk) yextreme(kk) tmin(ll) 1 1 1 0],[],[],[],[], ...
                            [1 1 1 1 1 1 -pi],[nX nY nT nX nY nT pi],[],options);
                    end
                end
                
                iminmse = find(mses == min(min(mses)));
                [kk,ll] = ind2sub(size(mses),iminmse);
                offargmin = squeeze(offargmins(kk,:,ll));
                
%                 onObjFn = objFnFn(maxZ);
%                 offObjFn = objFnFn(minZ);
%                 superObjFn = @(x) onObjFn(x(1:7)) + offObjFn(x([1 2 8:12]));
%                 argmin = fmincon(superObjFn, ...
%                     [xextreme(1) yextreme(1) tmax(1) 1 1 1 0 tmin(1) 1 1 1 0],[],[],[],[], ...
%                     [1 1 1 1 1 1 -pi 1 1 1 1 -pi],[nX nY nT nX nY nT pi nT nX nY nT pi],[],options);
                
%                 objFn = @(data,coeffFn) @(x) mean(mean((data(:,:,round(x(6)))-coeffFn(data)*gauss2d(X-x(4),Y-x(5),x(1),x(2),x(3),true)).^2));
%                 onObjFn = objFn(Z,@(x) max(max(max(x))));
%                 offObjFn = objFn(Z,@(x) min(min(min(x))));
%                 superObjFn = @(x) onObjFn(x(1:6)) + offObjFn(x([7:9 4 5 12])); % + abs(diff(x([6 12]))); % + ((x(4) - x(10)).^2 + (x(5) - x(11)).^2);
%                 
%                 optimalParams = fmincon(superObjFn,[1 1 0 xextreme(1) yextreme(1) textreme(1) 1 1 0 xextreme(1) yextreme(1) textreme(1)], ...
%                     [],[],[],[],[1 1 -pi 1 1 1 1 1 -pi 1 1 1],[nX nY pi nX nY nT nX nY pi nX nY nT],[],options);
                
%                 for kk = 1:numel(tmax)
%                     optimalParams = fmincon(objFn(Z,@(x) max(max(max(x)))),[1 1 0 xmax(kk) ymax(kk) tmax(kk)],[],[],[],[],[1 1 -pi 1 1 1],[nX nY pi nX nY nT],[],options);
                    
                    sdXOn = onargmin(4);
                    sdYOn = onargmin(5);
                    sdTOn = offargmin(6);
                    thetaOn = onargmin(7);
                    centreXOn = onargmin(1);
                    centreYOn = onargmin(2);
                    centreTOn = onargmin(3);
                    
                    if centreTOn > nT
                        centreTOn = nT;
                    end
                    
                    onRF = maxZ*gauss3d(X-centreXOn,Y-centreYOn,T-centreTOn,sdXOn,sdYOn,sdTOn,thetaOn,true);
                    
                    sdXOff = offargmin(4);
                    sdYOff = offargmin(5);
                    sdTOff = offargmin(6);
                    thetaOff = offargmin(7);
                    centreXOff = offargmin(1);
                    centreYOff = offargmin(2);
                    centreTOff = offargmin(3);
                    
                    if centreTOff > nT
                        centreTOff = nT;
                    end
                    
                    offRF = minZ*gauss3d(X-centreXOff,Y-centreYOff,T-centreTOff,sdXOff,sdYOff,sdTOff,thetaOff,true);
                    
                    RF = onRF + offRF;
                    
                    clf(fig);
                    subplot(2,2,1);
%                     figure;
                    hold on;
                    mesh(squeeze(Z(:,:,round(centreTOn))),'EdgeColor','k','FaceColor','none');
                    surf(RF(:,:,round(centreTOn)));
                    view([-37.5 30]);
%                     close gcf;
%                 end
                
%                 for kk = 1:numel(tmin)
%                     optimalParams = fmincon(objFn(Z,@(x) min(min(min(x)))),[1 1 0 xmin(kk) ymin(kk) tmin(kk)],[],[],[],[],[1 1 -pi 1 1 1],[nX nY pi nX nY nT],[],options);
                    
                    subplot(2,2,2);
%                     figure;
                    hold on;
                    mesh(squeeze(Z(:,:,round(centreTOff))),'EdgeColor','k','FaceColor','none');
                    surf(RF(:,:,round(centreTOff)));
                    view([-37.5 -30]);
                    
                    inOnRF = ((X-centreXOn)/(2*sdXOn)).^2 + ((Y-centreYOn)/(2*sdYOn)).^2 <= 1;
                    inOffRF = ((X-centreXOff)/(2*sdXOff)).^2 + ((Y-centreYOff)/(2*sdYOff)).^2 <= 1;
                    
%                     subplot(2,1,2);
%                     plot(squeeze(Z(round(mean([centreYOn centreYOff])),round(mean([centreXOn centreXOff])),:)));

% TODO : why doesn't this work?
%                     vid = figure;
%                     set(vid,'Position',[100 100 900 400]);                 
%                     
%                     for tt = 1:nT
%                         subplot(1,2,1);
%                         z = squeeze(Z(:,:,tt));
% %                         z(~(inOnRF | inOffRF)) = NaN;
%                         surf(z);
%                         caxis([minZ maxZ]);
%                         zlim([minZ maxZ]);
%                         colormap(gray);
%                         view(2);
%                         xlim([1 nX]);
%                         ylim([1 nY]);
%                         subplot(1,2,2);
%                         surf(squeeze(RF(:,:,tt)));
%                         caxis([minZ maxZ]);
%                         zlim([minZ maxZ]);
%                         colormap(gray);
%                         view(2);
%                         xlim([1 nX]);
%                         ylim([1 nY]);
% %                         input('...');
%                         frames(tt) = getframe(vid); %#ok<AGROW>
%                     end
%                         
%                     videoWriter = VideoWriter([saveFile '.avi']);
%                     videoWriter.FrameRate = 10;   
%                     open(videoWriter);
%                     writeVideo(videoWriter,frames);
%                     close(videoWriter);
%                     close(vid);
                    
                    temporalRF = zeros(nT,1);
                    
                    for tt = 1:nT
                        z = squeeze(Z(:,:,tt));
                        temporalRF(tt) = mean(z(inOnRF | inOffRF));
                    end
                    
                    subplot(2,1,2);
                    plot(1:nT,temporalRF);
                    
                    saveas(fig,[saveFile '.fig']);
                    saveas(fig,[saveFile '.png']);
%                 end
                
%                 if numel(tmax) > 1 || numel(tmin > 1)
%                     distances = sqrt( ...
%                         (repmat(ymax,1,numel(ymin))-repmat(ymin',numel(ymax),1)).^2 + ...
%                         (repmat(xmax,1,numel(xmin))-repmat(xmin',numel(xmax),1)).^2);
%                     
%                     iMinDistance = find(distances == min(min(distances)));
%                     
%                     if numel(iMinDistance) > 1
%                         error ffs;
%                     end
%                     
%                     [rMinDistance,cMinDistance] = ind2sub([tmax tmin],iMinDistance);
%                     
%                     ymax = ymax(rMinDistance);
%                     xmax = xmax(rMinDistance);
%                     tmax = tmax(rMinDistance);
%                     
%                     ymin = ymin(cMinDistance);
%                     xmin = xmin(cMinDistance);
%                     tmin = tmin(cMinDistance);
%                 end
%                     
%                 ymax = mean(ymax);
%                 xmax = mean(xmax);
%                 tmax = mean(tmax);
%                 ymin = mean(ymin);
%                 xmin = mean(xmin);
%                 tmin = mean(tmin);
%                 spatialMax = pixels(:,:,tmax);
%                 spatialMin = pixels(:,:,tmin);
%                 yy = reshape(Y,numel(Y),1);
%                 xx = reshape(X,numel(X),1);
%                 
%                 sdMax = std(reshape(spatialMax,numel(spatialMax),1));
%                 
%                 if sdMax == 0
%                     zmax = zeros(size(spatialMax));
%                 else
%                     zmax = (spatialMax - mean(reshape(spatialMax,numel(spatialMax),1)))/sdMax;
%                 end
%                 
%                 sdMin = std(reshape(spatialMin,numel(spatialMin),1));
%                     
%                 if sdMin == 0
%                     zmin = zeros(size(spatialMax));
%                 else
%                     zmin = (spatialMin - mean(reshape(spatialMin,numel(spatialMin),1)))/sdMin;
%                 end
                
%                 sigmays = 1:nY;
%                 sigmaxs = 1:nX;
%                 thetas = linspace(-pi,pi,100);
%                 
%                 for kk = 1 %:2
%                     Ny = numel(sigmays);
%                     Nx = numel(sigmaxs);
%                     mses = zeros(Ny,Nx,100);
%                     
%                     for ll = 1:Ny
%                         for mm = 1:Nx
%                             for oo = 1:numel(thetas)
%                                 tic;
%                                 offRF = gauss2d(X,Y,sigmaxs(mm),sigmays(ll),thetas(oo));
% %                                 offRF = mvnpdf([yy*sin(theta) xx*cos(theta)],[0 0],[sigmays(ll) sigmaxs(mm)]);
%                                 k = max(max(zmax))/max(max(offRF));
% %                                 offRF = reshape(offRF,size(zmax));
%                                 mses(ll,mm,oo) = mean(mean((zmax-k*offRF).^2));
%                                 toc;
%                             end
%                         end
%                     end
%                        
%                     index = find(mses == min(min(min(mses))));
%                     [y,x,z] = ind2sub(size(mses),index);
%                     
%                     sdy = sigmays(y);
%                     sdx = sigmaxs(x);
%                     phi = thetas(z);
%                     deltay = mean(diff(sigmays));
%                     deltax = mean(diff(sigmaxs));
%                     deltaz = mean(diff(thetas));
%                     
%                     sigmays = linspace(sdy - deltay,sdy + deltay,Ny);
%                     sigmaxs = linspace(sdx - deltax,sdx + deltax,Nx);
%                     thetas = linspace(phi - deltaz,phi+deltaz,100);
%                 end
% 
%                 clf;
% 
%                 objFn = @(data,coeffFn) @(x) mean(mean((data-coeffFn(data)*gauss2d(X-x(4),Y-x(5),x(1),x(2),x(3),true)).^2));
% %                 objFn = @(data,coeffFn,centreX,centreY) @(x) mean(mean((data-coeffFn(data)*gauss2d(X-centreX,Y-centreY,x(1),x(2),x(3),true)).^2));
%                 
%                 try
% %                     optimalParams = fmincon(objFn(zmax,@(x) max(max(x)),xmax,ymax),[1 1 0],[],[],[],[],[1 1 -pi],[nX nY pi],[],options);
%                 
%                     optimalParams = fmincon(objFn(zmax,@(x) max(max(x))),[1 1 0 xmax ymax],[],[],[],[],[1 1 -pi 1 1],[nX nY pi nX nY],[],options);
%                     sdXOn = optimalParams(1);
%                     sdYOn = optimalParams(2);
%                     thetaOn = optimalParams(3);
%                     centreXOn = optimalParams(4);
%                     centreYOn = optimalParams(5);
%                 catch err
%                     warning('On RF fitting failed for channel %s cluster index %d interation %d with error:\n',channelLabel,jj,nn); %#ok<WNTAG>
%                     logMatlabError(err);
%                     
%                     if nn > 1
%                         warning('Will use on RF from previous iteration'); %#ok<WNTAG>
%                     else
%                         warning('Aborting RF fitting for this cell...'); %#ok<WNTAG>
% %                         escape = true;
%                         break;
%                     end
%                 end
%                 
%                 onRF = max(max(zmax))*gauss2d(X-centreXOn,Y-centreYOn,sdXOn,sdYOn,thetaOn,true);
%                 subplot(2,2,1);
%                 hold on;
%                 mesh(zmax,'EdgeColor','k','FaceColor','none');
%                 surf(onRF);
% 
%                 try
% %                     optimalParams = fmincon(objFn(zmin,@(x) min(min(x)),xmin,ymin),[1 1 0],[],[],[],[],[1 1 -pi],[nX nY pi],[],options);
%                 
%                     optimalParams = fmincon(objFn(zmin,@(x) min(min(x))),[1 1 0 xmin ymin],[],[],[],[],[1 1 -pi 1 1],[nX nY pi nX nY],[],options);
%                     sdXOff = optimalParams(1);
%                     sdYOff = optimalParams(2);
%                     thetaOff = optimalParams(3);
%                     centreXOff = optimalParams(4);
%                     centreYOff = optimalParams(5);
%                 catch err
%                     warning('Off RF fitting failed for channel %s cluster index %d interation %d with error:\n',channelLabel,jj,nn); %#ok<WNTAG>
%                     logMatlabError(err);
%                     
%                     if nn > 1
%                         warning('Will use on RF from previous iteration'); %#ok<WNTAG>
%                     else
%                         warning('Aborting RF fitting'); %#ok<WNTAG>
% %                         escape = true;
%                         break;
%                     end
%                 end
%                 
%                 offRF = min(min(zmin))*gauss2d(X-centreXOff,Y-centreYOff,sdXOff,sdYOff,thetaOff,true);
%                 subplot(2,2,2);
%                 hold on;
%                 mesh(zmin,'EdgeColor','k','FaceColor','none');
%                 surf(offRF);
                
%                 offRF = offRF/max(max(offRF));
%                 offRF = reshape(k*offRF,size(zmax)); %#ok<NASGU>
                
%                 inOnRF = ((X-centreXOn)/(2*sdXOn)).^2 + ((Y-centreYOn)/(2*sdYOn)).^2 <= 1;
%                 inOffRF = ((X-centreXOff)/(2*sdXOff)).^2 + ((Y-centreYOff)/(2*sdYOff)).^2 <= 1;
% %                 tm = temporalRF;
%                 
%                 if ~any(any(inOnRF | inOffRF))
%                     % the receptive field is smaller than one pixel, so
%                     % it's probably not a real cell, in which case it
%                     % doesn't matter what we do
%                     temporalRF = ones(size(temporalRF));
%                     continue;
%                 end
% 
%                 for tt = 1:size(pixels,3)
%                     frame = pixels(:,:,tt);
%                     temporalRF(tt) = mean(frame(inOnRF | inOffRF));
%                 end
%                 
%                 temporalRF = temporalRF/sum(temporalRF);
%                 badPoints = temporalRF == 0 | isnan(temporalRF) | isinf(temporalRF);
%                 
%                 if any(badPoints)
%                     if all(badPoints)
%                         % I have no idea what this means, but it should
%                         % never happen for any sensible receptive field
%                         % profile
%                         temporalRF = ones(size(temporalRF));
%                     elseif any(badPoints(1:end-1) & badPoints(2:end))
%                         error('You''re in deep shit');
%                     else
%                         for badPoint = find(badPoints)
%                             if badPoint == 1 % extrapolate backwards
%                                 temporalRF(1) = 2*temporalRF(2)-temporalRF(3);
%                             elseif badPoint == numel(temporalRF) % extrapolate forwards
%                                 temporalRF(end) = 2*temporalRF(end-1)-temporalRF(end-2);
%                             else % interpolate
%                                 temporalRF(badPoint) = mean(temporalRF(badPoint+[-1 1]));
%                             end
%                         end
%                     end
%                 end
%                 
% %                 clf;
%                 subplot(2,1,2);
%                 plot(1:size(pixels,3),temporalRF);
%             end
            
%             if escape
%                 continue;
%             end
            
%             Zmax = Z(:,:,temporalRF == max(temporalRF));
%             Zmin = Z(:,:,temporalRF == min(temporalRF));
            Zmax = Z(:,:,round(centreTOn));
            Zmin = Z(:,:,round(centreTOff));
            Zmax = Zmax(inOnRF);
            Zmin = Zmin(inOffRF);
            pmax = normcdf(-sum(Zmax),0,numel(Zmax));
            pmin = normcdf(sum(Zmin),0,numel(Zmin));
%             pmax = 2*normcdf(-abs(sum(Zmax)),0,numel(Zmax));
%             pmin = 2*normcdf(-abs(sum(Zmin)),0,numel(Zmin));
            
            significant = ...
                   (~all(Zmax == mean(Zmax)) && pmax <= 0.05) ...
                || (~all(Zmin == mean(Zmin)) && pmin <= 0.05); %#ok<NASGU>
            
            save(saveFile,'RF','onRF','offRF','temporalRF','sdXOn','sdYOn','sdTOn','thetaOn','centreXOn','centreYOn','centreTOn','sdXOff','sdYOff','sdTOff','thetaOff','centreXOff','centreYOff','centreTOff','inOnRF','inOffRF','Zmax','Zmin','pmax','pmin','significant');
            fprintf('Done fitting RF for channel %s cluster #%d in %f seconds\n',channelLabel,jj,toc);
        end
    end
end

function fn = objFnFn(data,X,Y,T,k,useMSE)
    fn = @(x) objFn(x,data,X,Y,T,k,useMSE);
end

function fval = objFn(args,data,X,Y,T,k,useMSE)
    x0 = args(1);
    y0 = args(2);
    t0 = args(3);
    sdx = args(4);
    sdy = args(5);
    sdt = args(6);
    theta = args(7);
    
    if useMSE
        model = gauss3d(X-x0,Y-y0,T-t0,sdx,sdy,sdt,theta);

        fval = mean(mean(mean((data-k*model).^2)));
    else
        inRF = ((X-x0)/(2*sdx)).^2 + ((Y-y0)/(2*sdy)).^2 <= 1;
        
%         frame = data(:,:,round(t0))/sum(sum(inRF));
%         
%         fval = -abs(sum(frame(inRF)));
            frame = data(:,:,round(t0));
            frame = frame(inRF);
            fval = 2*normcdf(-abs(sum(frame)),0,numel(frame));
    end
end