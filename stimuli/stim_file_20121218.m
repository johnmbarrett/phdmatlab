repeats = 10;
speeds = 1:2;
angles = (0:11)*pi/6;
radius = 100;
centreX = 234+55;
centreY = 96+55;
maxT = 0;

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
            maxT = maxT + 2*radius*speed;
            stimulus = stimulus + 1;
            getExtraParams{stimulus} = struct( ...
                'startX',   startX,     ...
                'startY',   startY,     ...
                'finishX',  finishX,    ...
                'finishY',  finishY,    ...
                'length',   5*radius,   ...
                'maxT',     maxT,       ...
                'speed',    1/speed,    ...
                'width',    20);
            getPixels{stimulus} = @movingBar;
        end
    end
end

order = randperm(nStimuli);
getExtraParams2 = cell(size(getExtraParams));

for ii = 1:numel(getExtraParams)
    getExtraParams2{ii} = getExtraParams{order(ii)};
    getExtraParams2{ii}.maxT = getExtraParams{ii}.maxT;
end

stimulate(getPixels,[],[],[],60,'C:\John Barrett\JBWT0009\Stim 2.mat',false,false,getExtraParams2);