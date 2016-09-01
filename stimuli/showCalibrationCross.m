cross = calibrationCross(32,32,6);

myScreen = 0;
window = Screen('OpenWindow',myScreen,[],[1600 0 2240 480]);

tex = Screen('MakeTexture', window, cross);
Screen('DrawTexture', window, tex); 
Screen('Flip',window);

KbWait;

Screen('CloseAll');
clear Screen; 