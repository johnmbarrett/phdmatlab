maskTime = 60;
stimTime = 30;

width = 40; % 1000um / (25 px/um)
averageLuminance = 128;
contrasts = [1.0 0.4 0.2 0.1];
phases = (0:7)/8;

diffLuminances = round(averageLuminance*contrasts);
maxLuminances = min(averageLuminance+diffLuminances,255);
minLuminances = max(averageLuminance-diffLuminances,0);

disp(100*(maxLuminances-minLuminances)./(maxLuminances+minLuminances));

nContrasts = numel(contrasts);
nPhases = numel(phases);
nConditions = nWidths*nPhases;

nTrials = 150;

getPixels = cell(2*nConditions*nTrials+1,1);
textureRect = [120 0 640 480];

getPixels{1} = colourSequence(averageLuminance,textureRect,Inf,180);

for kk = 1:nTrials
    mask = colourSequence(averageLuminance,textureRect,Inf,maskTime);
    stims = cell(nConditions,1);

    for ii = 1:nContrasts
        for jj = 1:nPhases
            stims{sub2ind([nPhases nWidths],jj,ii)} = squareWave(maxLuminances(ii),width*2,Inf,pi/2,stimTime,minLuminances(ii),phases(jj));
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
markerRect = [0 0 textureRect(1) 480];

% !!! FILL IN SAVE FILE !!!
stimulate(getPixels,textureRect,[],markerRect,255,0,60,'D:\John B\Electrophysiology\Wild-Type Stim\JBWT0052\Stim 4.mat',false,false,false,cell(size(getPixels)));
    