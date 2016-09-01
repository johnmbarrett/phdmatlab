colour = 255;
colourss = {{colour 0*colour} {colour colour} {colour colour} {colour colour} {colour colour}};
rects = {[290 130 390 230;290 130 390 230] [290 130 340 230; 340 130 390 230] [340 130 390 230; 290 130 340 230] [290 130 390 180; 290 180 390 230] [290 180 390 230; 290 130 390 180]};
repeats = 15;
onTime = 2*60;
getPixels = cell(numel(colourss)*repeats,1);

for ii = 1:numel(colourss)
    f = colourSequence(colourss{ii},rects{ii},onTime/2,onTime);
    
    for jj = 1:repeats
        getPixels{repeats*(ii-1)+jj} = f;
    end
end

getPixels = getPixels(randperm(numel(getPixels)));

getPixels2 = cell(2*numel(getPixels),1);
getExtraParams = cell(size(getPixels2));

offTime = 5*60;
textureRect = [290 130 390 230];
black = colourSequence(0,textureRect,Inf,offTime);

for ii = 1:numel(getPixels)
    getPixels2{2*ii-1} = getPixels{ii};
    getPixels2{2*ii} = black;
    getExtraParams{2*ii-1} = struct([]);
    getExtraParams{2*ii} = struct([]);
end

trialTime = onTime+offTime;
totalTime = numel(getPixels) * trialTime;
disp(totalTime/(60*60));

% stimulate(getPixels2,textureRect,NaN,0,60,'D:\John B\Electrophysiology\Wild-Type Stim\JBWT0012\Stim 4.mat',false,false,getExtraParams);