% params
X = 334-(64-54)/2; % screen co-ords
Y = 247-(64-54)/2; % screen co-ords
repeats = 10;
colour = 255;
bgColour = 0;
screenRect = [];
textureRect = [0 0 640 480];
fps = 30;
squareTime = fps/2;
squareSize = 64/16;
stimulusMarkerRect = [0 0 120 480];
stimulusMarkerColour = 255;
whiteNoiseRect = [X Y X+64 Y+64];
whiteNoiseTime = 15*60*30; % frames
whiteNoiseParams = struct('isRandiDefined',true,'disableMarkerRect',true);
saveFile = 'D:\John B\Electrophysiology\Wild-Type Stim\JBWT0040\Stim 3.mat';

% stimulus generation code
quadrants = [0 0; 1 1; 0 1; 1 0];
mask = 3;
getPixels = cell(256,1);

for ii = 1:16
    for jj = 1:16
        index = sub2ind([16 16],jj,ii);
        x = 0;
        y = 0;

        for kk = 1:4
            bits = 2*(kk-1);
            quadIndex = bitshift(bitand(index-1,bitshift(mask,bits)),-bits);

            quad = quadrants(quadIndex+1,:);
            shift = 2^(4-kk);

            x = x + shift*quad(1);
            y = y + shift*quad(2);
        end

        getPixels{index} = colourSequence(colour,squareSize*[x y x+1 y+1]+[X Y X Y],Inf,squareTime);
    end
end

getPixels = repmat(getPixels,repeats,1);

getPixels{end+1} = whiteNoiseWithFrameMarker(whiteNoiseRect,stimulusMarkerRect,whiteNoiseTime,colour,squareSize,whiteNoiseParams);

disp(((size(getPixels,1)-1)*squareTime+whiteNoiseTime)/(60*fps));

stimulate(getPixels,textureRect,screenRect,stimulusMarkerRect,stimulusMarkerColour,bgColour,fps,saveFile,false,false,false,cell(size(getPixels)));
        