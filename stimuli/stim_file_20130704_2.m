X = 263;
Y = 233;
W = 112;
H = 112;
saveFile = 'D:\P25_4July13\P25_4July13_DA_WN.mat';

% !!! DON'T CHANGE ANYTHING BELOW THIS LINE !!!

% X = X-round((W-96)/2);
% Y = Y-round((H-96)/2);
screenRect = [1920 0 1920+640 480];
textureRect = [0 0 640 480];
fps = 30;
stimulusMarkerRect = [0 0 120 480];
stimulusMarkerColour = 255;
whiteNoiseRect = [X Y X+W Y+H];
whiteNoiseTime = 15; % minutes
whiteNoiseParams = struct('isRandiDefined',true);
colour = 255;
bgColour = 0;
pixelsPerSquare = 4;
getPixels = whiteNoiseWithFrameMarker(whiteNoiseRect,stimulusMarkerRect,whiteNoiseTime*60*fps,colour,pixelsPerSquare,whiteNoiseParams);
stimulate(getPixels,textureRect,screenRect,NaN,NaN,bgColour,fps,saveFile,false,false,false,struct([]));