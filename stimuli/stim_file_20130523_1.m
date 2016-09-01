% params
X = 265; % screen co-ords
Y = 235; % screen co-ords
width = 900; % microns
squareSizeUM = 107; % microns
pixelSizeUM = 26.79611650485436893203883495146;
xTiles = 3;
yTiles = 3;
repeats = 10;
colour = 255;
bgColour = 0;
whiteTime = 30;
blackTime = 30;
screenRect = [1920 0 1920+640 480];
textureRect = [0 0 640 480];
fps = 30;
stimulusMarkerRect = [0 0 120 480];
stimulusMarkerColour = 255;
saveFile = 'D:\P26_23May13\P26_23May13_LA_squares';

% stimulus generation code
squareSizePixels = ceil(squareSizeUM/pixelSizeUM);
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

disp((size(getPixels,1)*(whiteTime+blackTime))/(60*fps));

stimulate(getPixels2,textureRect,screenRect,stimulusMarkerRect,stimulusMarkerColour,bgColour,fps,saveFile,false,false,false,cell(size(getPixels2)));
        