maskTime = 60;
stimTime = 30;

width = 40; % 1000um / (25 px/um)
averageLuminance = 128;
luminances = [51 102 153 204];
phases = (0:7)/8;

nLuminances = numel(luminances);
nPhases = numel(phases);
nConditions = nLuminances*nPhases;

nTrials = 150;

getPixels = cell(2*nConditions*nTrials+1,1);
textureRect = [120 0 640 480];

getPixels{1} = colourSequence(averageLuminance,textureRect,Inf,180);

for kk = 1:nTrials
    mask = colourSequence(averageLuminance,textureRect,Inf,maskTime);
    stims = cell(nConditions,2);

    for ii = 1:nLuminances
        for jj = 1:nPhases
            row = sub2ind([nPhases nLuminances],jj,ii);
            stims{row,1} = squareWave(255,width*2,Inf,pi/2,stimTime,0,phases(jj));
            stims{row,2} = colourSequence(luminances(ii),textureRect,Inf,maskTime);
        end
    end

    stims = stims(randperm(nConditions),:);

    getPixels2 = cell(nConditions*2,1);

    for ii = 1:nConditions
        getPixels2{2*ii-1} = stims{ii,1};
        getPixels2{2*ii} = stims{ii,2};
    end
    
    getPixels(((kk-1)*nConditions*2+1:kk*nConditions*2)+1) = getPixels2;
end

disp((nConditions*nTrials*(maskTime+stimTime)+180)/(60*60));
markerRect = [0 0 textureRect(1) 480];

% !!! FILL IN SAVE FILE !!!
stimulate(getPixels,textureRect,[],markerRect,255,0,60,'D:\John B\Electrophysiology\Wild-Type Stim\JBWT0045\Stim 3.mat',false,false,false,cell(size(getPixels)));
    