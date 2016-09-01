function stimulate(getPixels, textureRect, screenRect, stimulusMarkerRect, stimulusMarkerColour, baseColour, fps, savefile, debug, leaveLastFrame, frameMarker, getExtraParams) 
    if nargin < 12 || isempty(getExtraParams) || ~(isa(getExtraParams, 'function_handle') || isstruct(getExtraParams) || iscell(getExtraParams));
        getExtraParams = @defaultGetExtraParams;
    end
    
    if nargin < 7 || ~isnumeric(fps) || ~isequal(size(fps),[1 1])
        fps = 60;
    end
    
    if nargin < 6 || isempty(baseColour)
        baseColour = [0 0 0];
    end
    
    if nargin < 5 || isempty(stimulusMarkerColour)
        stimulusMarkerColour = [255 255 255];
    end
    
    if numel(stimulusMarkerColour) ~= numel(baseColour)
        stimulusMarkerColour = 255*ones(size(baseColour))-baseColour;
    end
    
    if nargin < 4 || isempty(stimulusMarkerRect) || ~isnumeric(stimulusMarkerRect) || ~isreal(stimulusMarkerRect) || ~isequal(size(stimulusMarkerRect),[1 4]);
        stimulusMarkerRect = NaN;
    end

    if ~isnan(stimulusMarkerRect)
        stimulusMarkerWidth = diff(stimulusMarkerRect([1 3]));
        stimulusMarkerHeight = diff(stimulusMarkerRect([2 4]));
        stimulusMarkerPixels = zeros(stimulusMarkerWidth,stimulusMarkerHeight,numel(stimulusMarkerColour),2);
    
        for ii = 1:numel(stimulusMarkerColour)
            stimulusMarkerPixels(:,:,ii,1) = stimulusMarkerColour(ii)*ones(stimulusMarkerWidth,stimulusMarkerHeight);
            stimulusMarkerPixels(:,:,ii,2) = baseColour(ii)*ones(stimulusMarkerWidth,stimulusMarkerHeight);
        end
    end
    
    if nargin < 3 || isempty(screenRect) || ~isnumeric(screenRect) || ~isreal(screenRect) || ~isequal(size(screenRect),[1 4]);
        screenRect = [1920 0 1920+640 480]; % TODO : detect screen dimensions?
    end
    
    screenWidth = diff(screenRect([1 3]));
    screenHeight = diff(screenRect([2 4]));
    
    if nargin < 2 || isempty(textureRect) || ~isnumeric(textureRect) || ~isreal(textureRect) || ~isequal(size(textureRect),[1 4]);
        textureRect = [screenWidth/2-128 screenHeight/2-128 screenWidth/2+128 screenWidth/2+128];
    end
    
    textureWidth = diff(textureRect([1 3]));
    textureHeight = diff(textureRect([2 4]));
        
    if nargin < 1 || isempty(getPixels) || ~(isa(getPixels, 'function_handle') || iscell(getPixels));
        getPixels = @defaultGetPixels;
    end 
    
    nStimuli = size(getPixels,1);
    
    if iscell(getPixels) ~= iscell(getExtraParams) || numel(getPixels) ~= numel(getExtraParams)
        error 'Number of parameter sets does not match the number of stimuli';
    end
    
    stimulus = 1;
    if isstruct(getExtraParams)
        extraParams = getExtraParams;
    else
        extraParams = getExtraParams();
    end

    flipsPerFrame = floor(60/fps);
    texture = NaN;
    
    flipTimes = [];
    vbls = [];
    sots = [];
    ftss = [];
    misses = [];
    bposs = [];
    stimuli = [];
    
    versionStruct = ver('matlab');
    version = versionStruct.Version; %#ok<NASGU>
    
    if true || verLessThan('matlab','7.7')
        seed = rand('twister'); %#ok<NASGU,RAND>
    else
        seed = rng; %#ok<NASGU>
    end

    ListenChar(2);

    Screen('Preference', 'VisualDebugLevel', 1);
    myScreen = 0;
    window = Screen('OpenWindow',myScreen,baseColour,screenRect);
    
    if isnan(stimulusMarkerRect)
        stimulusMarkerTextures = nan;
    else
        stimulusMarkerTextures = zeros(2,1);
        
        for ii = 1:2
            stimulusMarkerTextures(ii) = Screen('MakeTexture', window, squeeze(stimulusMarkerPixels(:,:,:,ii)));
        end
    end
    
    % hold off starting the stimulation so I can have the shutter closed at
    % the beginning then open it just before I start stimulating
    KbWait;
    
    timeOffset = now*24*3600-GetSecs; %#ok<NASGU>
    
    escape = false;
    tt = 1;
    while tt < Inf
        tic;
        
        lastFrame = false;
        disableMarkerRect = false;
        if iscell(getPixels) && stimulus <= nStimuli
            pixelss = cell(1,size(getPixels,2));
            repeatFrames = zeros(1,size(getPixels,2));
            
            for ii = 1:size(getPixels,2);
                getPixel = getPixels{stimulus,ii};
                [pixels,extraParams,repeatFrame,markerIndex] = getPixel(tt,textureWidth,textureHeight,window,getExtraParams{stimulus});
                pixelss{ii} = pixels;
                lastFrame = lastFrame || any(any(any(isnan(pixels))));
                repeatFrames(ii) = repeatFrame;
                getExtraParams{stimulus,ii} = extraParams;
                disableMarkerRect = disableMarkerRect || (isfield(extraParams,'disableMarkerRect') && extraParams.disableMarkerRect);
            end
            
            repeatFrame = max(repeatFrames);
        else
            [pixels,extraParams,repeatFrame] = getPixels(tt,textureWidth,textureHeight,window,extraParams);
            pixelss = {pixels};
            lastFrame = any(any(isnan(pixels)));
            disableMarkerRect = isfield(extraParams,'disableMarkerRect') && extraParams.disableMarkerRect;
        end
        
        if lastFrame
            if ~iscell(getPixels) || stimulus == nStimuli
                break;
            else
                tt = 1;
                stimulus = stimulus + 1;
                
                continue;
            end
        end
        if ~isnan(repeatFrame)
            repeats = repeatFrame * flipsPerFrame;
        else
            repeats = flipsPerFrame;
        end
%         fprintf('Stimulus %d generation iteration: %d time: %f s\n',stimulus,tt,toc);
        
        tic;
        if ~isnan(texture) && texture > 0
            Screen('Close',texture);
            texture = NaN;
        end    

        for ii = 1:numel(pixelss)
            pixels = pixelss{ii};
            
            if ~isempty(pixels)
                texture = Screen('MakeTexture', window, pixels);

                % TODO : multiple textureRects for multiple simultaneous
                % stimuli
                Screen('DrawTexture', window, texture, [], textureRect); 
            end
        end
        
        if ~any(isnan(stimulusMarkerRect)) && ~disableMarkerRect
            if isnan(markerIndex)
                if frameMarker
                    markerIndex = ceil(tt/(repeats/flipsPerFrame));
    %                 fprintf('%d\t%d\t%d\n',tt,repeats,markerIndex);
                else
                    markerIndex = stimulus;
                end
            end

            Screen('DrawTexture', window, stimulusMarkerTextures(mod(markerIndex-1,2)+1), [], stimulusMarkerRect);
        end 
        
%         fprintf('Texture manipulation iteration: %d time: %f s\n',tt,toc);
        
        for jj = 1:repeats
            tic;
            [keyIsDown,secs,keyCode] = KbCheck(-1); %#ok<ASGLU>

            if keyIsDown
                if nargin > 6 && all(debug)
                    disp(KbName(keyCode));
                end

                if keyCode(KbName('esc'))
                    escape = true;
                    break;
                end
            end
%             fprintf('Keyboard poll iteration: %d time: %f s\n',jj,toc);

            [vbl, sot, fts, miss, bpos] = Screen('Flip',window,0,double(jj < repeats || (leaveLastFrame && stimulus == nStimuli)));
            flipTimes = [flipTimes; now]; %#ok<AGROW>
            
            if nargin > 8 && all(debug)
                fprintf('%f\t%f\t%f\t%f\t%f\n', vbl, sot, fts, miss, bpos); 
            end
            
            vbls = [vbls; vbl]; %#ok<AGROW>
            sots = [sots; sot]; %#ok<AGROW>
            ftss = [ftss; fts]; %#ok<AGROW>
            misses = [misses; miss]; %#ok<AGROW>
            bposs = [bposs; bpos]; %#ok<AGROW>
            stimuli = [stimuli; stimulus]; %#ok<AGROW>
        end
        
        if escape
            break;
        end

        if ~isnan(repeatFrame)
            tt = tt + repeatFrame;
        else
            tt = tt + 1;
        end
    end

    % end of experiment tidy up
    if ~any(isnan(stimulusMarkerRect))
        for ii = 1:2
            Screen('Close',stimulusMarkerTextures(ii));
        end
    end

    if ~isnan(texture) && texture > 0
        Screen('Close',texture);
    end

    Screen('Flip',window);
    
    KbWait;
    
    ListenChar(0);
    Screen('CloseAll');

    clear Screen; 
    
    if nargin > 7 && ischar(savefile) && ~isempty(savefile)
        save(savefile, 'getPixels', 'getExtraParams', 'vbls', 'sots', 'ftss', 'misses', 'bposs', 'stimuli', 'timeOffset', 'seed', 'extraParams', 'version', 'textureRect', 'screenRect', 'baseColour', 'fps', 'flipTimes');
    end
end

function extraParams = defaultGetExtraParams()
    extraParams = struct('isRandiDefined',~verLessThan('matlab','7.7'),'maxT',18000);
end