function uledSTA(retina,stimulus,detectionRecordings,staRecordings)
    if nargin < 1
        retina = 'A';
    end
    
    if nargin < 2
        stimulus = 'flashing squares';
%         stimulus = 'moving bars';
    end
    
    if nargin < 3
        detectionRecordings = [3 11];
    end
    
    if isnumeric(retina)
        retina = char('A'+retina-1);
    end
    
    %%

    sequence = load([stimulus ' sequence.mat']);

    if strcmp(stimulus,'flashing squares')
        timingVariable = 'duration';
        conditions = [sequence.conditions ones(size(sequence.conditions,1),1)];
        timingIndex = size(conditions,2);
        subFrames = 1;
        ctrlRecording = 6;
        drugRecording = 16;

        onPixelsFunction = @(X,Y,condition,sequence,varargin) ...
            X >= sequence.X(condition(:,1)) & ...
            X < sequence.X(condition(:,1))+sequence.width(condition(:,3)) & ...
            Y >= sequence.Y(condition(:,2)) & ...
            Y < sequence.Y(condition(:,2))+sequence.width(condition(:,3));

        skipFunction = @(condition) condition(:,3) == 1;
    elseif strcmp(stimulus,'moving bars')
        timingVariable = 'period';
        timingIndex = 2;
        subFrames = 16-sequence.width;
        ctrlRecording = 7;
        drugRecording = 17;

        angles = [4 3 2 1 5 6 7 0]*pi/4;

        onPixelsFunction = @(X,Y,condition,sequence,subFrame) ...
            X*cos(angles(condition(:,1))) + Y*sin(angles(condition(:,1))) >= subFrame-1 & ...
            X*cos(angles(condition(:,1))) + Y*sin(angles(condition(:,1))) < subFrame+width-1;
    end
    
    if nargin > 4
        ctrlRecording = staRecordings(1);
        drugRecording = staRecordings(2);
    end

    recordings = [ctrlRecording drugRecording];

    nStimuli = numel(sequence.conditionOrder);
    Ts = sequence.(timingVariable);
    [Y,X] = ndgrid(0:15,0:15);

    %%

    load(['Retina ' retina '_stimulusTimes.mat']);

    stimulusTimes = [onsetss{recordings}]; %#ok<USENS>
    
    if ~all(isfinite(stimulusTimes(:)))
        warning('Stimulus times missing or corrupt, aborting...\n');
        return;
    end

    load(['Retina ' retina '_spiketimestamps.mat']);

    valid = cells(:,2) > 0; %#ok<NODEF>
    spiketimestamps = spiketimestamps(valid); %#ok<NODEF>
    cells = cells(valid,:);
    nCells = size(cells,1);

    %%

    stas = zeros(16,16,2,nCells);
    ns = zeros(16,16);

    for ii = 1:nStimuli
        condition = conditions(sequence.conditionOrder(ii),:);

        if skipFunction(condition)
            continue;
        end

        tic;
        T = Ts(condition(timingIndex))/1000;

        for jj = 1:subFrames
            pixels = zeros(16,16);
            pixels(onPixelsFunction(X,Y,condition,sequence,jj)) = 1;
            ns = ns + pixels;

            for kk = 1:nCells
                spikeTimes = spiketimestamps{kk};

                for ll = 1:2
                    t = stimulusTimes(ii,ll);
                    nSpikes = sum(spikeTimes > t & spikeTimes < t+T);
                    stas(:,:,ll,kk) = stas(:,:,ll,kk) + pixels*nSpikes/T;
                end
            end
        end
        toc;
    end

    stas = stas./repmat(ns,[1 1 2 nCells]);

    %%

    load(['Retina ' retina '_uled_responsive_cells.mat']);
    responsiveCells = intersect(responsive{detectionRecordings}); %#ok<USENS>
%     nResponses = numel(responsiveCells);

    %%

    xvals = [reshape(X,256,1) reshape(Y,256,1)]+1;
    gaussFun = @(params,xvals) params(1)*gauss2d(xvals(:,1)-params(2),xvals(:,2)-params(3),params(4),params(5),params(6),true)+params(7);

    paramss = zeros(7,2,nCells);
    residuals = zeros(nCells,2);

    options = optimset('lsqcurvefit');
    
    %%

    for ii = 1:nCells
        tic;
        for jj = 1:2
            sta = reshape(stas(:,:,jj,ii),256,1);

            if std(sta) == 0
                paramss(:,jj,ii) = [0 8.5 8.5 Inf Inf 0 0]';
                continue;
            end

            sta = sta./sum(sta);
            [maxS,maxI] = max(sta);
            [maxY,maxX] = ind2sub([16 16],maxI(1));
            [paramss(:,jj,ii),residuals(ii,jj)] = lsqcurvefit(gaussFun,[maxS maxX maxY 1 1 0 mean(sta)],xvals,sta,[0 1 1 0 0 0 0],[1 16 16 8 8 2*pi maxS],options);
        end
        toc;
    end
    
    %%
    
    amplitude = squeeze(paramss(1,:,:))'; %#ok<NASGU>
    centreX = squeeze(paramss(2,:,:))'; %#ok<NASGU>
    centreY = squeeze(paramss(3,:,:))'; %#ok<NASGU>
    sdX = squeeze(paramss(4,:,:))'; %#ok<NASGU>
    sdY = squeeze(paramss(5,:,:))'; %#ok<NASGU>
    angle = squeeze(paramss(6,:,:))'; %#ok<NASGU>
    baseline = squeeze(paramss(7,:,:))'; %#ok<NASGU>

    %%

    titles = {'Control' 'Drug'};

    figure
    set(gcf,'Position',[0 0 1200 1200],'Renderer','zbuffer');

    padsurf = @(X,varargin) surf([X X(:,end); X(end,:) X(end,end)],varargin{:});

    for ii = 1:numel(responsiveCells)
        responseIndex = responsiveCells(ii);
        
        sta = stas(:,:,:,responseIndex);
        cc = [min(sta(:)) max(sta(:))];

        for jj = 1:2
            subplot(2,2,jj);
            padsurf(sta(:,:,jj));
            shading interp;
            
            if cc(1) < cc(2)
                caxis(cc);
            end
                
            title(titles{jj});
            view(2);
            xlim([1 16]);
            ylim([1 16]);
            set(gca,'XTick',[],'YTick',[]);

            if jj == 1
                ylabel('Data');
            end

            subplot(2,2,2+jj);
            padsurf(reshape(gaussFun(paramss(:,jj,responseIndex),xvals),16,16));
            shading interp;
            view(2);
            xlim([1 16]);
            ylim([1 16]);
            set(gca,'XTick',[],'YTick',[]);

            if jj == 1
                ylabel('Gaussian Fit');
            end
        end

        channel = cells(responseIndex,1);
        cluster = cells(responseIndex,2);
        suptitle(sprintf('Channel %d cluster %d',channel,cluster));

        figFile = sprintf('Retina %c_sta_%s_channel_%d_cluster_%d',retina,strrep(stimulus,' ','_'),channel,cluster);
        saveas(gcf,figFile,'fig');
        saveas(gcf,figFile,'png');
    end

    close(gcf);

    %%

    % stamean = zeros(nCells,2);
    % stastd = zeros(nCells,2);
    % 
    % for ii = 1:nCells
    %     for jj = 1:2
    %         x0 = paramss(2,jj,ii);
    %         y0 = paramss(3,jj,ii);
    %         sdx = paramss(4,jj,ii);
    %         sdy = paramss(5,jj,ii);
    %         phi = paramss(6,jj,ii);
    %         
    %         Xcan = (X-x0)*cos(phi) + (Y-y0)*sin(phi);
    %         Ycan = (Y-y0)*cos(phi) - (X-x0)*sin(phi);
    %         inRF = (Xcan/(2*sdx)).^2+(Ycan/(2*sdy)).^2 <= 1;
    %         
    %         sta = stas(:,:,jj,ii);
    %         stamean(ii,jj) = mean(sta(inRF));
    %         stastd(ii,jj) = std(sta(~inRF));
    %     end
    % end
    % 
    % stasnr = stamean./stastd;

    %%

    save(['Retina ' retina '_uled_' strrep(stimulus,' ','_') 'sta.mat'],'stas','amplitude','centreX','centreY','sdX','sdY','angle','baseline','residuals');
end