% params
X = 233; % screen co-ords
Y = 145; % screen co-ords
width = 1000; % microns
squareSizeUM = 100; % microns
xTiles = 2;
yTiles = 2;
repeats = 10;
colour = 255;
bgColour = 0;
whiteTime = 30;
blackTime = 30;
screenRect = [];
textureRect = [0 0 640 480];
fps = 30;
stimulusMarkerRect = [520 0 640 480];
stimulusMarkerColour = 255;
whiteNoiseRect = [243-8 155-8 243-8+96 155-8+96];
whiteNoiseTime = 15; % minutes
whiteNoiseParams = struct('isRandiDefined',false);
saveFile = 'C:\John Barrett\JBWT0028\Stim 1.mat';

% stimulus generation code
whiteNoiseParams.bgColour = bgColour;
whiteNoiseParams.disableMarkerRect = true;
whiteNoiseTime = whiteNoiseTime*fps*60;
squareSizePixels = ceil(squareSizeUM/25);
nSquares = ceil(width/squareSizeUM);
W = squareSizePixels*nSquares;
nTiles = xTiles*yTiles;
tileCorners = repmat([X Y],nTiles,1)+[(repmat(W*(0:xTiles-1)',yTiles,1)) (kron(W*(0:yTiles-1)',ones(xTiles,1)))];

getPixels = cell(repeats*nSquares^2,nTiles);
allSquares = [repmat((1:nSquares)',nSquares,1) kron((1:nSquares)',ones(nSquares,1))];
randSquares = allSquares(randperm(nSquares^2),:);
xs = randSquares(:,1);
ys = randSquares(:,2);

for jj = 1:nSquares^2
    for ii = 1:nTiles
        x = tileCorners(ii,1);
        y = tileCorners(ii,2);
        squareRect = [x y x y] + squareSizePixels*[xs(jj)-1 ys(jj)-1 xs(jj) ys(jj)];
        index = (1:repeats)+repeats*(jj-1);
        getPixels(index,ii) = repmat({colourSequence(colour,squareRect,Inf,whiteTime)},repeats,1);
    end
end
        
getPixels2 = cell(size(getPixels).*[2 1]);
getPixels2(1:2:end-1,:) = getPixels;

for ii = 2:2:size(getPixels2,1)
    for jj = 1:nTiles
        x = tileCorners(jj,1);
        y = tileCorners(jj,2);
        getPixels2{ii,jj} = colourSequence(bgColour,[x y x+W y+W],Inf,blackTime);
    end
end

getPixels2(end+1,:) = {whiteNoiseWithFrameMarker(whiteNoiseRect,stimulusMarkerRect,whiteNoiseTime,colour,2,whiteNoiseParams) @emptyStimulus @emptyStimulus @emptyStimulus};

disp((size(getPixels,1)*(whiteTime+blackTime)+whiteNoiseTime)/(60*fps));

stimulate(getPixels2,textureRect,screenRect,stimulusMarkerRect,stimulusMarkerColour,bgColour,fps,saveFile,false,false,false,cell(size(getPixels2)));
        