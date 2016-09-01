function stimulate2(getPixels, textureRect, screenRect, stimulusMarkerRect, baseColour, fps, savefile, debug, leaveLastFrame, getExtraParams) 
    if nargin < 9 || isempty(getExtraParams) || ~(isa(getExtraParams, 'function_handle') || isstruct(getExtraParams) || iscell(getExtraParams));
        getExtraParams = @defaultGetExtraParams;
    end
    
    if nargin < 5 || ~isnumeric(fps) || ~isequal(size(fps),[1 1])
        fps = 60;
    end
    
    if nargin < 4 || isempty(baseColour)
        baseColour = [0 0 0];
    end
    
    if nargin < 3 || isempty(screenRect) || ~isnumeric(screenRect) || ~isreal(screenRect) || ~isequal(size(screenRect),[1 4]);
        screenRect = [1680 0 1680+640 480]; % TODO : detect screen dimensions?
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
    
    if iscell(getPixels) ~= iscell(getExtraParams) || numel(getPixels) ~= numel(getExtraParams)
        error 'Number of parameter sets does not match the number of stimuli';
    end
    
    stimulus = 1;
    if isstruct(getExtraParams)
        extraParams = getExtraParams;
    else
        extraParams = getExtraParams();
    end

    repeats = floor(60/fps);
    texture = NaN;
    ii = 1;
    tt = 1;
    
    flipTimes = [];
    vbls = [];
    sots = [];
    ftss = [];
    misses = [];
    bposs = [];
    stimuli = [];
    
    versionStruct = ver('matlab');
    version = versionStruct.Version; %#ok<NASGU>
    
    if verLessThan('matlab','7.7')
        seed = rand('twister'); %#ok<NASGU,RAND>
    else
        seed = rng; %#ok<NASGU>
    end

    ListenChar(2);

    Screen('Preference', 'VisualDebugLevel', 1);
    myScreen = 0;
    window = Screen('OpenWindow',myScreen,baseColour,screenRect);
    
    % hold off starting the stimulation so I can have the shutter closed at
    % the beginning then open it just before I start stimulating
    KbWait;
    
    timeOffset = now*24*3600-GetSecs; %#ok<NASGU>
    
    while ii < Inf
        [keyIsDown,secs,keyCode] = KbCheck(-1); %#ok<ASGLU>
    
        if keyIsDown
            if nargin > 6 && all(debug)
                disp(KbName(keyCode));
            end
        
            if keyCode(KbName('esc'))
                break;
            end
        end
        
        if iscell(getPixels) && stimulus <= numel(getPixels)
            getPixel = getPixels{stimulus};
            [pixels,extraParams] = getPixel(tt,textureWidth,textureHeight,window,getExtraParams{stimulus});
            getExtraParams{stimulus} = extraParams;
        else
            [pixels,extraParams] = getPixels(tt,textureWidth,textureHeight,window,extraParams);
        end
        
        if any(isnan(pixels))
            if ~iscell(getPixels) || stimulus == numel(getPixels)
                break;
            else
                tt = 1;
                stimulus = stimulus + 1;
                continue;
            end
        end

        for jj = 1:repeats
            if ~isnan(texture) && texture > 0
                Screen('Close',texture);
            end

            if ~isempty(pixels)
                texture = Screen('MakeTexture', window, pixels);

                Screen('DrawTexture', window, texture, [], textureRect); 
            end

            [vbl, sot, fts, miss, bpos] = Screen('Flip',window);
            flipTimes = [flipTimes; now];
            
            if nargin > 6 && all(debug)
                fprintf('%f\t%f\t%f\t%f\t%f\n', vbl, sot, fts, miss, bpos); 
            end
            
            vbls = [vbls; vbl]; %#ok<AGROW>
            sots = [sots; sot]; %#ok<AGROW>
            ftss = [ftss; fts]; %#ok<AGROW>
            misses = [misses; miss]; %#ok<AGROW>
            bposs = [bposs; bpos]; %#ok<AGROW>
            stimuli = [stimuli; stimulus]; %#ok<AGROW>
        end

        ii = ii + 1;
        tt = tt + 1;
    end

    if ~leaveLastFrame && ~isnan(texture)
        % blank the screen at the end of the experiment rather than leaving
        % the last stimulus frame on screen
        Screen('Close',texture);
        Screen('Flip',window);
    end
    
    KbWait;
    
    ListenChar(0);
    Screen('CloseAll');

    clear Screen; 
    
    if nargin > 5 && ischar(savefile) && ~isempty(savefile)
        save(savefile, 'getPixels', 'getExtraParams', 'vbls', 'sots', 'ftss', 'misses', 'bposs', 'stimuli', 'timeOffset', 'seed', 'extraParams', 'version', 'textureRect', 'screenRect', 'baseColour', 'fps', 'flipTimes');
    end
end

function extraParams = defaultGetExtraParams()
    extraParams = struct('isRandiDefined',~verLessThan('matlab','7.7'),'maxT',18000);
end