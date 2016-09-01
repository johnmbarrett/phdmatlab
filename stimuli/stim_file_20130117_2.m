whiteTime = 1*60;
blackTime = 3*60;
radii = (50:50:600)/(2*25); % (diameter in microns)/(2*microns per pixel) = radius in pixels
centres = [315 155; 365 155; 315 205; 365 205];
repeats = 10;
getPixels = cell(numel(radii),size(centres,1),repeats);

for hh = 1:numel(radii)
    for ii = 1:size(centres,1)
        f = spot(centres(ii,:),radii(hh),255,whiteTime);

        for jj = 1:repeats
            getPixels{hh,ii,jj} = f;
        end
    end
end

getPixels = getPixels(randperm(numel(getPixels)));

getPixels2 = cell(2*numel(getPixels),1);
getExtraParams = cell(size(getPixels2));

black = colourSequence(0,[0 0 640 480],Inf,blackTime);

for ii = 1:numel(getPixels)
    getPixels2{2*ii-1} = getPixels{ii};
    getPixels2{2*ii} = black;
    getExtraParams{2*ii-1} = struct([]);
    getExtraParams{2*ii} = struct([]);
end

trialTime = blackTime+whiteTime;
totalTime = numel(getPixels) * trialTime;
disp(totalTime/(60*60));

textureRect = [290 130 390 230];
stimulate(getPixels2,textureRect,NaN,0,60,'D:\John B\Electrophysiology\Wild-Type Stim\JBWT0012\Stim 5.mat',false,false,getExtraParams);
% stimulate(getPixels2,textureRect,NaN,0,60,NaN,false,false,getExtraParams);