[filename, filepath] = uigetfile('*.mj2');

video = VideoReader([filepath filename]);
N = video.NumberOfFrames;
%%
disp('GET TO DA CHOPPA!!!');
% Screen('Preference', 'SkipSyncTests', 2);
myScreen = 0;
window = Screen('OpenWindow',myScreen,[],[1600 0 2400 600  ]);

fps = 60;
spf = 1/fps;
repeats = floor(60/fps);
tnext = 0;
wold = NaN;
for ii = 1:N
    frame = read(video,ii);

    for jj = 1:repeats
        if ~isnan(wold)
            Screen('Close',wold);
        end

        wnew = Screen('MakeTexture', window, squeeze(frame));

        Screen('DrawTexture', window, wnew); 
    
        [vbl, sot, fts] = Screen('Flip',window);%,tnext);
        fprintf('%f %f %f\n', vbl, sot, fts); 
%          tnext = tcurr + spf;
    end

    wold = wnew; 
end

KbWait;
Screen('CloseAll');

clear Screen; 