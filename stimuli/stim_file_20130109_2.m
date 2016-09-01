periods = 2.^(0:2:6);
repeats = 5;
maxT = max(periods)*5;
textureRects = [0 0 279 480; 0 0 640 213; 279 0 640 480; 0 213 640 480];
getPixels = cell(size(textureRects,1),numel(periods),repeats);

for hh = 1:size(textureRects,1)
    textureRect = textureRects(hh,:);
    
    for ii = 1:numel(periods)
        f = colourSequence({255,0},textureRect,periods(ii),maxT);

        for jj = 1:repeats
            getPixels{hh,ii,jj} = f;
        end
    end
end

getPixels = getPixels(randperm(numel(getPixels)));

getPixels2 = cell(2*numel(getPixels),1);
getExtraParams = cell(size(getPixels2));

trialTime = 15*60;
blackTime = trialTime - maxT;
black = colourSequence(0,textureRect,Inf,blackTime);

for ii = 1:numel(getPixels)
    getPixels2{2*ii-1} = getPixels{ii};
    getPixels2{2*ii} = black;
    getExtraParams{2*ii-1} = struct([]);
    getExtraParams{2*ii} = struct([]);
end

totalTime = numel(getPixels) * trialTime;
disp(totalTime/(60*60));

stimulate(getPixels2,textureRect,NaN,0,60,'C:\John Barrett\JBWT0011\Stim 2.mat',false,false,getExtraParams);