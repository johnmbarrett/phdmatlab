% cross = checkBoard(64,64,1,1);
cross = zeros(64,64,3);
cross(:,:,[1 3]) = 255*ones(64,64,2);

myScreen = 0;
window = Screen('OpenWindow',myScreen,[],[1600 0 2240 480]);

tex = Screen('MakeTexture', window, cross);
Screen('DrawTexture', window, tex); 
Screen('Flip',window);

KbWait;

Screen('CloseAll');
clear Screen; 