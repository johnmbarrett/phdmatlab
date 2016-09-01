periods = 2.^(4:6);
repeats = 10;
maxT = max(periods)*5;
textureRect = [0 0 640 480];
getPixels = cell(numel(periods)*repeats,1);

for ii = 1:numel(periods)
    f = colourSequence({255,0},textureRect,periods(ii),maxT);
    
    for jj = 1:repeats
        getPixels{repeats*(ii-1)+jj} = f;
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

stimulate(getPixels2,textureRect,NaN,0,60,'D:\John B\Electrophys\Testing\Photodiode 01\Stim 6.mat',false,false,getExtraParams);