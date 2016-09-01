function [concData,timeData,voltData,infoData,ns] = aggregateStimData(experiments,drugs,stimOffsets,middleVolts)
    if nargin < 3
        stimOffsets = zeros(size(experiments));
    end
    
    nExperiments = numel(experiments);
    nDrugs = numel(drugs);
    
    ns = nan(nExperiments,nDrugs,6,3);
    
    if nargin < 4
        middleVolts = 8*ones(nExperiments,nDrugs);
    end
    
    infoData = nan(nExperiments,nDrugs,2,2);
    sheets = {'Count' 'Latency'};
    
    for ii = 1:nExperiments
        expDir = sprintf('JBOG%04d',experiments(ii));
        
        for jj = 1:nDrugs
            infoFile = sprintf('%s\\%s_%s_info.xlsx',expDir,expDir,drugs{jj});
            
            if ~exist(infoFile,'file')
                continue;
            end
            
            excelData = importdata(infoFile);
            
            for kk = 1:2
                sheetData = excelData.data.(sheets{kk});
                
                nCols = size(sheetData,2);
                nReps = (nCols-1)/2;
                
                sheetData = permute(reshape(sheetData(:,2:end),[2 nReps 2]),[3 2 1]);
                
                infoData(ii,jj,:,kk) = squeeze(median(diff(sheetData)));
            end
        end
    end
    
    data = cell(nExperiments,nDrugs);
    
    concStims = [3 7:9 12 16 19 23:25 28 32];
    voltStims = [1:5 10:14 17:21 26:30];
    allStims = union(concStims,voltStims);
    nStims = numel(allStims);
    
    sigmoid = @(b,x) 1./(1+exp(-(x-b(1))/b(2)));
    beta0 = [0 1];
    pws = [5;10;25;50;75;100];
    
    % can't use a for loop here due to the edge case where one experiment
    % has the same drug twice and MATLAB for loops don't work like C ones
    ii = 1;
    while ii <= nExperiments
        expDir = sprintf('JBOG%04d',experiments(ii));
        responsiveCellss = cell(nDrugs,1);
        
        for jj = 1:nDrugs
            responsiveCellFile = sprintf('%s\\%s_%s_responsive_cells.mat',expDir,expDir,drugs{jj});
            
            if exist(responsiveCellFile,'file')
                responsiveCells = load(responsiveCellFile,'responsiveCells');
                responsiveCellss{jj} = responsiveCells.responsiveCells;
                nCells = size(responsiveCells,1);
                data{ii,jj} = zeros(nCells,6,5,6,4);
            end
        end
        
        recDirs = dir(sprintf('%s\\Stim*',expDir));
        recDirs = {recDirs(vertcat(recDirs.isdir)).name}';
        
        seenDrugs = {};
        
        for jj = 1:nStims
            stim = allStims(jj);
            stimName = sprintf('Stim %d ',stim+stimOffsets(ii));
            
            recIndex = strncmpi(stimName,recDirs,numel(stimName));
            
            if ~any(recIndex)
                continue;
            end
            
            recDir = recDirs{recIndex};
            
            tokens = regexp(recDir,[stimName '([0-9]+) ([a-zA-Z0-9]+) ([0-9]+)V'],'tokens','once');
            conc = str2double(tokens{1});
            drug = tokens{2};
            volt = str2double(tokens{3});
            
            if mod(stim,16) == 1 && ismember(drug,seenDrugs)
                data = [data(1:ii,:); cell(1,nDrugs); data((ii+1):end,:)];
                experiments = experiments([1:ii ii (ii+1):end]);
                stimOffsets = stimOffsets([1:ii ii (ii+1):end]);
                middleVolts = middleVolts([1:ii ii (ii+1):end],:);
                ns = [ns(1:ii,:,:,:); nan(1,nDrugs,6,3); ns((ii+1):end,:,:,:)];
                ii = ii + 1;
                nExperiments = nExperiments+1;
            else
                seenDrugs{end+1} = drug; %#ok<AGROW>
            end
            
            dataFile = sprintf('%s\\%s\\%s_uled_square_responses_newmethod.mat',expDir,recDir,recDir);
            
            if ~exist(dataFile,'file')
                continue;
            end
            
            load(dataFile);
            n = size(responsiveCells,1);
            
            drugIndex = find(strcmpi(drug,drugs));
            
            if isempty(drugIndex)
                continue;
            end
            
            if isempty(responsiveCellss{drugIndex})
                drugSuffix = ceil(stim/16);
                responsiveCellFile = sprintf('%s\\%s_%s_%d_responsive_cells.mat',expDir,expDir,drugs{drugIndex},drugSuffix);
                
                if ~exist(responsiveCellFile,'file')
                    continue;
                end
                
                responsiveCells = load(responsiveCellFile,'responsiveCells');
                responsiveCells = responsiveCells.responsiveCells;
                
                data{ii,drugIndex} = zeros(nCells,6,5,6,4);
            else
                responsiveCells = responsiveCellss{drugIndex};
            end
            
            nCells = size(responsiveCells,1);
            
            if conc == 0 
                concIndex = 1+5*(mod(stim,16) == 0);
            else
                concIndex = 2+log(conc/10)/log(2);
            end
            
            voltIndex = volt-5;
            
            if volt == middleVolts(ii,drugIndex)
                ns(ii,drugIndex,concIndex,1) = n;
            end
            
            if ismember(concIndex,[1 5])
                ns(ii,drugIndex,voltIndex,2+(concIndex>1)) = n;
            end
            
            responseIndices = [];
            cellIndices = [];
            
            for kk = 1:nCells
                channel = responsiveCells(kk,1);
                cluster = responsiveCells(kk,2);
                
                index = find(ismember([channels clusters],[channel cluster],'rows'));
                
                if ~isempty(index)
                    responseIndices(end+1) = kk; %#ok<AGROW>
                    cellIndices(end+1) = index; %#ok<AGROW>
                end
            end
            
            fiddle = @(f,x) permute(f(x),[3 2 1]);
                
            nSpikes = allNSpikes(:,:,cellIndices); %#ok<NODEF>
            mSpikes = fiddle(@median,nSpikes);
            
            data{ii,drugIndex}(responseIndices,:,voltIndex,concIndex,1) = mSpikes;
            
            dSpikes = mSpikes-repmat(meds(cellIndices),1,6);
            
            data{ii,drugIndex}(responseIndices,:,voltIndex,concIndex,2) = dSpikes;
            
            sigmas(sigmas == 0) = sqrt(1/1200); %#ok<AGROW>
            zSpikes = (fiddle(@mean,nSpikes)-repmat(mus(cellIndices),1,6))./repmat(sigmas(cellIndices),1,6);
%             zSpikes = fiddle(@mean,nSpikes)./fiddle(@std,nSpikes);
            
            data{ii,drugIndex}(responseIndices,:,voltIndex,concIndex,3) = zSpikes;
            
            pSpikes = allPSpikes(:,cellIndices); %#ok<NODEF>
            thresholds = zeros(numel(cellIndices),1);
            
            for kk = 1:size(pSpikes,3)
                bp = nlinfit(pws,pSpikes(:,kk),sigmoid,beta0);
                thresholds(kk) = min(200,max(0,bp(1))); % I need to learn about censoring
            end
            
            data{ii,drugIndex}(responseIndices,6,voltIndex,concIndex,4) = thresholds;
            
%             data{ii,drugIndex}(:,:,voltIndex,concIndex,2) = snr(:,:)'; %#ok<NODEF>
%             data{ii,drugIndex}(:,6,voltIndex,concIndex,3) = thresholds;
        end
        
        ii = ii + 1;
    end
    
    concData = nan(nExperiments,6,nDrugs,4);
    timeData = nan(nExperiments,6,2,nDrugs,4);
    voltData = nan(nExperiments,5,2,nDrugs,4);
    
    for ii = 1:nExperiments    
        for jj = 1:nDrugs
            X = data{ii,jj};
            
            if isempty(X)
                continue;
            end
            
            voltIndex = middleVolts(ii,jj)-5;
            
            for kk = 1:4
                valid = X(:,6,voltIndex,1,kk) ~= 0;
                C = permute(X(valid,6,voltIndex,:,kk),[1 4 2 3 5]);
                C = C./repmat(C(:,1),1,6);
                C(isnan(C)) = 1; % Inf/Inf = 1 here
                concData(ii,:,jj,kk) = nanmedian(C,1);
                
                if kk < 4
                    valid = X(:,6,3,1,kk) ~= 0;
                    T = permute(X(valid,:,3,[1 5],kk),[1 2 4 3 5]);
                    T = T./repmat(T(:,6,1),[1 6 2]);
                    timeData(ii,:,:,jj,kk) = nanmedian(T,1);
                end
                
                valid = X(:,6,3,1,kk) ~= 0;
                V = permute(X(valid,6,:,[1 5],kk),[1 3 4 2 5]);
                V = V./repmat(V(:,3,1),[1 5 2]);
                V(isnan(V)) = 1;
                voltData(ii,:,:,jj,kk) = nanmedian(V,1);
            end
        end
    end
end