W = 800;
H = 600;
screenRect = [1920 0 1920+W H];
luminances = [0 15:16:255];
n = numel(luminances);

% Screen('Preference', 'SkipSyncTests', 2);
Screen('Preference', 'VisualDebugLevel', 1);
myScreen = 0;
window = Screen('OpenWindow',myScreen,[0 0 0],screenRect);

ListenChar(2);


ii = 0;
oldTex = NaN;
doneOnce = false;
while true
    [keyIsDown,blah,keyCode] = KbCheck(-1);

    if ~doneOnce || keyIsDown
        doneOnce = true;

        if keyCode(KbName('esc'))
            break;
        end

        ii = ii + 1;
        
        index = mod(ii-1,n)+1;
    
        pixels = luminances(index)*ones(W,H);

        if ~isnan(oldTex)
            Screen('Close',oldTex);
        end

        newTex = Screen('MakeTexture', window, pixels);

        Screen('DrawTexture', window, newTex, [], [0 0 W H]); 

        Screen('Flip',window);
        
        pause(1);
         
        oldTex = newTex;
    end
end

ListenChar(0);
Screen('CloseAll');

clear Screen; 