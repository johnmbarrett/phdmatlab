% params
X = 290-(64-60)/2; % screen co-ords
Y = 210-(64-60)/2; % screen co-ords
repeats = 10;
colour = 255;
bgColour = 0;
screenRect = [];
textureRect = [120 0 640 480];
fps = 30;
pauseTime = 8*fps;
flashTime = 2*fps;
squareTime = fps/2;
squareSize = 64/16;
stimulusMarkerRect = [0 0 120 480];
stimulusMarkerColour = 255;
whiteNoiseRect = [X Y X+64 Y+64];
whiteNoiseTime = 20*60*fps; % frames
whiteNoiseParams = struct('isRandiDefined',true,'disableMarkerRect',true);
saveFile = 'D:\John B\Electrophysiology\Wild-Type Stim\JBWT0054\Stim 1.mat';

% stimulus generation code
getPixels1 = cell(4,1);

% using measured monitor intensities, 136 gives a brightness about halfway
% between 0 and 255
getPixels1{1} = colourSequence(136,textureRect,pauseTime,pauseTime);
getPixels1{2} = colourSequence({255; 136; 0; 136},textureRect,repmat([flashTime; pauseTime],2,1),repeats*(flashTime+pauseTime)*2,true,1);
getPixels1{3} = colourSequence(0,textureRect,pauseTime,pauseTime/2,true,1);
getPixels1{4} = colourSequence(0,textureRect,pauseTime,pauseTime/2);

quadrants = [0 0; 1 1; 0 1; 1 0];
mask = 3;
getPixels2 = cell(256,1);

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

        getPixels2{index} = colourSequence(colour,squareSize*[x y x+1 y+1]+[X Y X Y],Inf,squareTime);
    end
end

getPixels2 = repmat(getPixels2,repeats,1);

getPixels = [getPixels1; getPixels2];

getPixels{end+1} = whiteNoiseWithFrameMarker(whiteNoiseRect,stimulusMarkerRect,whiteNoiseTime,colour,squareSize,whiteNoiseParams);

disp((2*pauseTime+repeats*(pauseTime+flashTime)*2+(size(getPixels2,1)-1)*squareTime+whiteNoiseTime)/(60*fps));

stimulate(getPixels,[0 0 640 480],screenRect,stimulusMarkerRect,stimulusMarkerColour,bgColour,fps,saveFile,false,false,false,cell(size(getPixels)));
        