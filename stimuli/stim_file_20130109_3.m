trialTime = 5*60;
whiteTime = 60;
blackTime = trialTime-whiteTime;
radii = (50:50:800)/(2*25); % (diameter in microns)/(2*microns per pixel) = radius in pixels
centres = [230 165; 230 215; 280 165; 280 215];
repeats = 5;
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

totalTime = numel(getPixels) * trialTime;
disp(totalTime/(60*60));

stimulate(getPixels2,textureRect,NaN,0,60,'C:\John Barrett\JBWT0011\Stim 3.mat',false,false,getExtraParams);
% stimulate(getPixels2,textureRect,NaN,0,60,NaN,false,false,getExtraParams);