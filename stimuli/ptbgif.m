[filename, filepath] = uigetfile('*.gif');

m = imread([filepath filename],'gif','frames','all');
%%
disp('GET TO DA CHOPPA!!!');
% Screen('Preference', 'SkipSyncTests', 2);
myScreen = 0;
window = Screen('OpenWindow',myScreen);% ,[],[100 100 256 256]);
white = WhiteIndex(window); % pixel value for white
black = BlackIndex(window); % pixel value for black
gray = (white+black)/2;
inc = white-gray;
Screen(window, 'FillRect', gray);

tnext = 0;
fps = 10;
spf = 1/fps;
wold = NaN;
for ii = 1:size(m,4)
%     w = Screen(window, 'MakeTexture', gray+inc*squeeze(m(:,:,1,ii))/255);
    if ~isnan(wold)
        Screen('Close',wold);
    end
    
    wnew = Screen('MakeTexture', window,  squeeze(m(:,:,1,ii)));

    Screen('DrawTexture', window, wnew);
    tcurr = Screen('Flip',window,tnext);
    tnext = tcurr + spf;
     
    wold = wnew;
end

KbWait;
Screen('CloseAll');

clear Screen; 