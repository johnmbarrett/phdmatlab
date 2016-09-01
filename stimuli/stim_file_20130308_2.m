colour = 255;
repeats = 7;
speeds = 1:2;
angles = (0:11)*pi/6;
radius = 80;
centreX = 240;
centreY = 160;

nStimuli = repeats*length(speeds)*length(angles);
getPixels = cell(nStimuli,1);
getExtraParams = cell(nStimuli,1);

stimulus = 0;
for ii = 1:repeats
    for jj = 1:2
        speed = speeds(jj);
        
        for kk = 1:12
            angle = angles(kk);
            
            startX = centreX-radius*cos(angle);
            startY = centreY-radius*sin(angle);
            finishX = centreX+radius*cos(angle);
            finishY = centreY+radius*sin(angle);
            stimulus = stimulus + 1;
            getExtraParams{stimulus} = struct( ...
                'colour',   colour,         ...
                'maxT',     2*radius*speed, ...
                'startX',   startX,         ...
                'startY',   startY,         ...
                'finishX',  finishX,        ...
                'finishY',  finishY,        ...
                'length',   5*radius,       ...
                'speed',    1/speed,        ...
                'width',    20);
            getPixels{stimulus} = @movingBar;
        end
    end
end

getExtraParams2 = getExtraParams(randperm(numel(getExtraParams)));

disp(numel(getExtraParams2)*2*radius*speed/(60*60));

stimulate(getPixels,[],[],[530 0 640 480],colour,0,60,'C:\John Barrett\JBWT0021\Stim 9.mat',false,false,false,getExtraParams2);
% stimulate(getPixels,[],[],NaN,NaN,0,60,'NaN',false,false,false,getExtraParams2);