colour = 255;
colourss = {{colour 0*colour} {colour colour} {colour colour} {colour colour} {colour colour}};

X = 182;
Y = 127;
W = 92;
H = 92;
rects = { ...
    [X Y X+W Y+H; X Y X+W Y+H] ...
    [X Y X+W/2 Y+H; X+W/2 Y X+W Y+H] ...
    [X+W/2 Y X+W Y+H; X Y X+W/2 Y+H] ...
    [X Y X+W Y+H/2; X Y+H/2 X+W Y+H] ...
    [X Y+H/2 X+W Y+H; X Y X+W Y+H/2]};

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
textureRect = [X Y X+W Y+H];
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

stimulate(getPixels2,textureRect,[],[560 0 640 480],colour,zeros(size(colour)),60,'C:\John Barrett\JBWT0016\Stim 1.mat',false,false,getExtraParams);
% stimulate(getPixels2,textureRect,[],[560 0 640 480],colour,zeros(size(colour)),60,NaN,false,false,getExtraParams);