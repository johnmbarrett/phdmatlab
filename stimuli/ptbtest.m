myScreen = 0;
window = Screen(myScreen, 'OpenWindow');
white = WhiteIndex(window); % pixel value for white
black = BlackIndex(window); % pixel value for black
gray = (white+black)/2;
inc = white-gray;
Screen(window, 'FillRect', gray);
[x,y] = meshgrid(-100:100, -100:100);
g = exp(-((x/50).^2)-((y/50).^2));

w = zeros([size(g) 100]);

for i = 1:100
    n = randi(2,size(g))-1;
    m = conv2(n,g,'same');
    w(:,:,i) = Screen(window, 'MakeTexture', gray+inc*m);
end

for i = 1:100
    Screen('DrawTexture', window, w(i));
    Screen(window,'Flip'); 
end

KbWait;
Screen('CloseAll');

clear Screen; 