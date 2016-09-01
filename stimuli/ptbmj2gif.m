[gifname, gifpath] = uigetfile('*.gif');
[filename, filepath] = uigetfile('*.mj2');

video = VideoReader([filepath filename]);
N = video.NumberOfFrames;

gif = imread([gifpath gifname],'frames','all');
%%
disp('GET TO DA CHOPPA!!!');
% Screen('Preference', 'SkipSyncTests', 2);
myScreen = 0;
window = Screen('OpenWindow',myScreen);% ,[],[100 100 256 256]);

% movie = Screen('OpenMovie',window,[filepath filename]);
% Screen('PlayMovie',movie,1);

wold1 = NaN;
wold2 = NaN;
for ii = 1:N
    frame1 = read(video,ii);
    frame2 = gif(:,:,1,ii);
    
    if ~isnan(wold1) && ~isnan(wold2)
        Screen('Close',wold1);
        Screen('Close',wold2);
    end
    
    wnew1 = Screen('MakeTexture', window, squeeze(frame1));
    wnew2 = Screen('MakeTexture', window, squeeze(frame2));

    Screen('DrawTexture', window, wnew1, [], [0 0 256 256]);
    Screen('DrawTexture', window, wnew2, [], [256 0 512 256]);
    tcurr = Screen('Flip',window); %%v ,tnext);
%     tnext = tcurr + spf;
     
    wold1 = wnew1;
    wold2 = wnew2;
end

KbWait;
Screen('CloseAll');

clear Screen; 