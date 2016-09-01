maskTime = 60;
stimTime = 30;

widths = 2.^(3:6);
phases = (0:7)/8;

nWidths = numel(widths);
nPhases = numel(phases);
nConditions = nWidths*nPhases;

nTrials = 50;

getPixels = cell(2*nConditions*nTrials+1,1);
textureRect = [0 0 520 480];

getPixels{1} = colourSequence(64,textureRect,Inf,180);

for kk = 1:nTrials
    mask = colourSequence(64,textureRect,Inf,maskTime);
    stims = cell(nConditions,1);

    for ii = 1:nWidths
        for jj = 1:nPhases
            stims{sub2ind([nPhases nWidths],jj,ii)} = squareWave(128,widths(ii)*2,Inf,pi/2,stimTime,0,phases(jj));
        end
    end

    stims = stims(randperm(numel(stims)));

    getPixels2 = cell(nConditions*2,1);

    for ii = 1:nConditions
        getPixels2{2*ii-1} = stims{ii};
        getPixels2{2*ii} = mask;
    end
    
    getPixels(((kk-1)*nConditions*2+1:kk*nConditions*2)+1) = getPixels2;
end

disp((nConditions*nTrials*(maskTime+stimTime)+180)/(60*60));
markerRect = [textureRect(3) 0 640 480];

% !!! FILL IN SAVE FILE !!!
stimulate(getPixels,textureRect,[],markerRect,255,0,60,'C:\John Barrett\JBWT0027\Stim 6.mat',false,false,false,cell(size(getPixels)));
    