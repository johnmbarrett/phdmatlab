nStimuli = 1+4*4*4+8*5*5;
getPixels = cell(1,nStimuli);
getExtraParams = cell(1,nStimuli);

getPixels{1} = @defaultGetPixels;
maxT = 15*60*60;
getExtraParams{1} = struct('isRandiDefined',false,'maxT',maxT);

textureRect = [277 103 105 105];
midX = textureRect(1) + textureRect(3)/2;
midY = textureRect(2) + textureRect(4)/2;
startXOffset = [1 -1 0 0];
startYOffset = [0 0 1 -1];
finishXOffset = [-1 1 0 0];
finishYOffset = [0 0 -1 1];

stimulus = 1;
for hh = 1:4
    speed = 2^(hh-1);
    
    for ii = 1:4
        width = 2^(ii-1);

        for jj = 1:4
            stimulus = stimulus + 1;
            getPixels{stimulus} = @movingBar;
            maxT = maxT + floor(textureRect(3)/speed);

            getExtraParams{stimulus} = struct(...
                'startX',   midX + startXOffset(jj)*textureRect(3)/2, ...
                'startY',   midX + startYOffset(jj)*textureRect(4)/2, ...
                'finishX',  midX + finishXOffset(jj)*textureRect(3)/2, ...
                'finishY',  midX + finishYOffset(jj)*textureRect(4)/2, ...
                'width',    width, ...
                'speed',    speed, ...
                'maxT',     maxT);
        end
    end
end

for hh = 1:8
    angle = (hh-1)*pi/4;
    
    for ii = 1:5
        spatialPeriod = 2^(ii+2);
        
        for jj = 1:5
            temporalPeriod = 2^(ii+2);
            
            stimulus = stimulus+1;
            getPixels{stimulus} = @grating;
            
            maxT = maxT + 10*temporalPeriod;
            
            getExtraParams{stimulus} = struct(...
                'angle',            angle, ...
                'spatialPeriod',    spatialPeriod, ...
                'temporalPeriod',   temporalPeriod, ...
                'maxT',             maxT);
        end
    end
end

disp(maxT/(60*60));