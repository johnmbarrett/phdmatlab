function epileptogenesis(framerate,colour,screenRect) 
    if nargin < 3
        screenRect = [1680 0 1680+640 480];
    end

    if nargin < 2
        colour = 255;
    end
    
    if nargin < 1
        framerate = 60;
    end
    
    pixelss = zeros(480,640,numel(colour),2);
    
    for ii = 1:numel(colour)
        pixelss(:,:,ii,2) = colour(ii)*ones(480,640);
    end

    Screen('Preference', 'VisualDebugLevel', 1);
    myScreen = 0;
    window = Screen('OpenWindow',myScreen,0,screenRect);
    
    escape = false;
    repeats = ceil(60/framerate);
    texture = NaN;
    ii = 0;
    
    % hold off starting the stimulation so I can have the shutter closed at
    % the beginning then open it just before I start stimulating
    KbWait;
    
    while ii < Inf
        pixels = squeeze(pixelss(:,:,:,mod(ii,2)+1));

        for jj = 1:repeats
            [keyIsDown,secs,keyCode] = KbCheck(-1); %#ok<ASGLU>

            if keyIsDown
                if keyCode(KbName('esc'))
                    escape = true;
                    break;
                end
            end
            
            if ~isnan(texture) && texture > 0
                Screen('Close',texture);
            end

            texture = Screen('MakeTexture', window, pixels);

            Screen('DrawTexture', window, texture,[],[0 0 640 480]);
        
            Screen('Flip',window);
        end
        
        ii = ii + 1;
        
        if escape
            break;
        end
    end

    KbWait;
    
    ListenChar(0);
    Screen('CloseAll');

    clear Screen;
end