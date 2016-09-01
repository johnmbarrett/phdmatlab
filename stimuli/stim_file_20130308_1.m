%%
colour = 255;
whiteTime = 60;
blackTime = 60;
radii = (400:200:800)'/(2*25); % (diameter in microns)/(2*microns per pixel) = radius in pixels
X = 189;
Y = 112;
W = 104;
X = X+floor((W-96)/2);
Y = Y+floor((W-96)/2);
W = 96;
[centreX,centreY] = meshgrid(1:5,1:5);
centreXY = W*[reshape(centreX,25,1) reshape(centreY,25,1)]/6;
centres = repmat([X Y],25,1) + centreXY;
repeats = 10;
N = numel(radii)*size(centres,1)*repeats;
centreIdx = repmat((1:25)',numel(radii)*repeats,1);
radiusIdx = repmat(kron((1:3)',ones(size(centres,1),1)),repeats,1);

%%
remainingIndices = 1:N;
invalidIndices = {};
badDraws = 0;
permutedIndices = nan(N,1);

for ii = 1:N
    invalid = unique([invalidIndices{:}]);
    valid = setdiff(remainingIndices,invalid);
    
    if ~isempty(valid)
        nextIndex = valid(randi(numel(valid)));
    else
        nextIndex = remainingIndices(randi(numel(remainingIndices)));
        badDraws = badDraws + 1;
    end
    
    permutedIndices(ii) = nextIndex;
    
    remainingIndices = setdiff(remainingIndices,nextIndex);
    
    if isempty(remainingIndices)
        break;
    end
    
    XY = centres(centreIdx(nextIndex),:);
    radius = radii(radiusIdx(nextIndex));
    
    newlyInvalid = remainingIndices(sqrt(sum((repmat(XY,numel(remainingIndices),1)-centres(centreIdx(remainingIndices),:)).^2,2)) <= radius+radii(radiusIdx(remainingIndices)));
    
    invalidIndices{end+1} = newlyInvalid; %#ok<SAGROW>
    
    if numel(invalidIndices) > 1
        invalidIndices = invalidIndices(2:end);
    end
end
 
disp(badDraws);

%%
% while ~validPermutation
%     tic;
%     idx = randperm(N);
%     validPermutation = true;
%     
%     for ii = 1:N-5
%         x1 = centreXY(centreIdx(idx(ii)),1);
%         y1 = centreXY(centreIdx(idx(ii)),2);
%         r1 = radii(radiusIdx(idx(ii)));
%         
%         for jj = 1:5
%             x2 = centreXY(centreIdx(idx(ii+jj)),1);
%             y2 = centreXY(centreIdx(idx(ii+jj)),2);
%             r2 = radii(radiusIdx(idx(ii+jj)));
%             
%             if sqrt((x1-x2).^2 + (y1-y2).^2) <= (r1+r2)
%                 validPermutation = false;
%                 break;
%             end
%         end
%         
%         if ~validPermutation
%             break;
%         end
%     end
%     toc;
% end
% 
% return;
%%
getPixels = cell(N,1);
getExtraParams = cell(N,1);
for ii = 1:N
    getPixels{ii} = spot(centres(centreIdx(permutedIndices(ii)),:),radii(radiusIdx(permutedIndices(ii))),colour,whiteTime);
    getExtraParams{ii} = struct([]);
end

%%
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

%%
textureRect = [X Y X+W Y+W];
stimulate(getPixels2,textureRect,NaN,[530 0 640 480],colour,zeros(size(colour)),60,'C:\John Barrett\JBWT0021\Stim 1.mat',false,false,false,getExtraParams);
% stimulate(getPixels2,textureRect,NaN,[560 0 640 480],colour,zeros(size(colour)),60,NaN,false,false,falsegetExtraParams);
% stimulate(getPixels2,textureRect,[],[560 0 640 480],colour,zeros(size(colour)),60,NaN,false,false,false,getExtraParams);