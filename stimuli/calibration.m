function values = calibration(colour,screenRect,fps)
    if nargin < 3
        fps = 6;
    end
    
    if nargin < 2 || isempty(screenRect) || any(isnan(screenRect))
        screenRect = [1920 0 1920+640 480];
    end
    
    width = 60;
    height = 60;
    x = (640-width)/2;
    y = (480-height)/2;
    inc = 1;
    luminance = 128;
    contrast = 100;

    % Screen('Preference', 'SkipSyncTests', 2);
    Screen('Preference', 'VisualDebugLevel', 1);
    myScreen = 0;
    window = Screen('OpenWindow',myScreen,[0 0 0],screenRect);

    ListenChar(2);


    ii = 0;
    oldTex = NaN;
    while true
        [keyIsDown,blah,keyCode] = KbCheck(-1);

        if keyIsDown
    %         disp(KbName(keyCode));

            if keyCode(KbName('esc'))
                break;
            end

            if keyCode(KbName('left'))
                x = x - inc;
            end

            if keyCode(KbName('right'))
                x = x + inc;
            end

            if keyCode(KbName('up'))
                y = y - inc;
            end

            if keyCode(KbName('down'))
                y = y + inc;
            end

            if keyCode(KbName('=+'))
                width = width + inc;
                height = height + inc;
            end

            if keyCode(KbName('-_'))
                width = max(1,width - inc);
                height = max(1,height - inc);
            end

            if keyCode(KbName('-'));
                inc = max(1,inc - 1);
            end

            if keyCode(KbName('+'))
                inc = inc + 1;
            end

            if keyCode(KbName('q'))
                contrast = min(100,contrast + 1); 
            end

            if keyCode(KbName('a'))
                contrast = max(0,contrast - 1); 
            end

            if keyCode(KbName('w'))
                luminance = min(255,luminance + 1); 
            end

            if keyCode(KbName('s'))
                luminance = max(0,luminance - 1); 
            end
        end

        if mod(ii,ceil(60/fps)) == 0
            I = 255*ones(width,height);
            Z = zeros(width,height);
            noise = min(I,max(Z,luminance + luminance*contrast*(2*randint(width,height)-1)/100));
            
            if nargin > 0
                n = numel(colour);
                pixels = zeros([size(noise) n]);
                for ii = 1:n
                    pixels(:,:,ii) = colour(ii)*noise/255;
                end
            else
                pixels = noise;
            end
        end

        if ~isnan(oldTex)
            Screen('Close',oldTex);
        end

        newTex = Screen('MakeTexture', window, pixels);

        Screen('DrawTexture', window, newTex, [], [x y x+width y+height]); 

        [vbl, sot, fts] = Screen('Flip',window);
    %     fprintf('%f %f %f\n', vbl, sot, fts); 

        oldTex = newTex; 
        ii = ii+1;
    end

    ListenChar(0);
    Screen('CloseAll');

    clear Screen; 

    values = [x y width height luminance contrast];
end