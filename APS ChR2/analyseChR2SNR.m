% TODO : double-dipping

load('channelNames.mat');
load('ChR2_cells.mat');
load('EventNo.mat');
load('spiketimestamps.mat');

%%

[nCells,nRecs] = size(allNSpikes);
snr = zeros(nCells,nRecs);

for ii = 1:nCells
    tic;
    for jj = 1:nRecs
        sample = bootstrap(allNSpikes{ii,jj}{2},30);
        
        if isequal(sample,zeros(size(sample)))
            sample(1) = 1;
        end
        
        sigma = std(sample/0.1);
        snr(ii,jj) = peakFR(ii,jj,1)/sigma;
    end
    toc;
end
        