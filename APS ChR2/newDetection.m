if ~exist('detectionFile','var')
    detectionFile = 'ChR2_cells_ff.mat';
end

load(detectionFile);

%%

if ~exist('spontStartIndex','var')
    spontStartIndex = 16;
end

%%

nTrials = size(nSpikesss{1}{1},1);

nRecs = size(nSpikesss,2);
allBetterUnits = cell(1,4);

for ii = 1:nRecs
    allBetterUnits{ii} = find(cellfun(@(c) max(cellfun(@(A) mean(sum(A(:,1:5),2),1),c)),nSpikesss(:,ii)) >= 1);
end

%%

beforeSamples = cell(nRecs,2);
afterMaxFR = cell(nRecs,2);
pChange = cell(nRecs,2);
allBestUnits = cell(nRecs,2);

for ii = 1:nRecs
    for jj = 1:2
        nCells = numel(allBetterUnits{ii});
        disp(nCells);
        
        beforeSamples{ii,jj} = cell(nCells,1);
        afterMaxFR{ii,jj} = zeros(nCells,1);
        pChange{ii,jj} = zeros(nCells,1);
        
        for kk = 1:nCells
            tic;
            cellIndex = allBetterUnits{ii}(kk);
            beforeSamples{ii,jj}{kk} = bootstrap(nSpikesss{cellIndex,ii}{3-jj}(1:2:end,spontStartIndex:end),ceil(nTrials/2));
            afterMaxFR{ii,jj}(kk) = max(mean(nSpikesss{cellIndex,ii}{jj}(1:2:end,1:5),1));
            pLeft = sum(beforeSamples{ii,jj}{kk} <= afterMaxFR{ii,jj}(kk))/numel(beforeSamples{ii,jj}{kk});
            pRight = sum(beforeSamples{ii,jj}{kk} >= afterMaxFR{ii,jj}(kk))/numel(beforeSamples{ii,jj}{kk});
            pChange{ii,jj}(kk) = 2*min(pLeft,pRight);
            toc;
        end
        
        h = pChange{ii,jj} <= 0.05/numel(allBetterUnits{ii});
        allBestUnits{ii,jj} = goodUnits(allBetterUnits{ii}(h));
    end
end

%%

if ~exist('saveFile','var')
    saveFile = 'ChR2_cells_new.mat';
end

save(saveFile,'-v7.3','afterMaxFR','allBestUnits','allBetterUnits','beforeSamples','goodUnits','nSpikesss','pChange','samples')