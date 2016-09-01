board = checkBoard(640,480,32,32);
crossLum = calibrationCross(32,32,6);
crossRGB = zeros([size(crossLum) 3]);
crossRGB(:,:,1) = 255-crossLum;
crossRGB(:,:,3) = 255-crossLum;

Screen('Preference', 'VisualDebugLevel', 1);
myScreen = 0;
window = Screen('OpenWindow',myScreen,[],[1920 0 1920+640 480]); %[1680 0 2320 480]);

boardTex = Screen('MakeTexture', window, board);
crossTex = Screen('MakeTexture', window, crossRGB);
Screen('DrawTexture', window, boardTex); 
Screen('DrawTexture', window, crossTex); 
Screen('Flip',window);

KbWait;

Screen('CloseAll');
clear Screen; 