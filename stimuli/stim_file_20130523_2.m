screenRect = [1920 0 1920+640 480];
textureRect = [0 0 640 480];
fps = 30;
stimulusMarkerRect = [0 0 120 480];
stimulusMarkerColour = 255;
whiteNoiseRect = [294 221 294+96 221+96];
whiteNoiseTime = 15; % minutes
whiteNoiseParams = struct('isRandiDefined',true);
colour = 255;
bgColour = 0;
pixelsPerSquare = 4;
saveFile = 'C:\John Barrett\JBWT0029\Stim 1.mat';
getPixels = whiteNoiseWithFrameMarker(whiteNoiseRect,stimulusMarkerRect,whiteNoiseTime*60*fps,colour,pixelsPerSquare,whiteNoiseParams);
stimulate(getPixels,textureRect,screenRect,NaN,NaN,bgColour,fps,saveFile,false,false,false,struct([]));