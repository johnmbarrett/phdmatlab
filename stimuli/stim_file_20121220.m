repeats = 4;
speeds = 1:2;
angles = (0:11)*pi/6;
radius = 250;
centreX = 194+111;
centreY = 67+111;

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

order = randperm(nStimuli);
getExtraParams2 = cell(size(getExtraParams));

for ii = 1:numel(getExtraParams)
    getExtraParams2{ii} = getExtraParams{order(ii)};
end

stimulate(getPixels,[],[],[],60,'C:\John Barrett\JBWT0010\Stim 2.mat',false,false,getExtraParams2);