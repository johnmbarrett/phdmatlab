colour = 255;
whiteTime = 1*60;
blackTime = 3*60;
radii = (100:50:700)/(2*25); % (diameter in microns)/(2*microns per pixel) = radius in pixels
X = 182;
Y = 127;
W = 93;
centres = repmat([X Y],4,1) + W*[1 1; 2 1; 1 2; 2 2]/3 ;
repeats = 10;
getPixels = cell(numel(radii),size(centres,1),repeats);

for hh = 1:numel(radii)
    for ii = 1:size(centres,1)
        f = spot(centres(ii,:),radii(hh),colour,whiteTime);

        for jj = 1:repeats
            getPixels{hh,ii,jj} = f;
        end
    end
end

getPixels = getPixels(randperm(numel(getPixels)));

getPixels2 = cell(2*numel(getPixels),1);
getExtraParams = cell(size(getPixels2));

black = colourSequence(0,[0 0 560 480],Inf,blackTime);

for ii = 1:numel(getPixels)
    getPixels2{2*ii-1} = getPixels{ii};
    getPixels2{2*ii} = black;
    getExtraParams{2*ii-1} = struct([]);
    getExtraParams{2*ii} = struct([]);
end

trialTime = blackTime+whiteTime;
totalTime = numel(getPixels) * trialTime;
disp(totalTime/(60*60));

textureRect = [X Y X+W Y+W];
stimulate(getPixels2,textureRect,NaN,[560 0 640 480],colour,zeros(size(colour)),60,'C:\John Barrett\JBWT0016\Stim 2.mat',false,false,getExtraParams);
% stimulate(getPixels2,textureRect,NaN,[560 0 640 480],colour,zeros(size(colour)),60,NaN,false,false,getExtraParams);