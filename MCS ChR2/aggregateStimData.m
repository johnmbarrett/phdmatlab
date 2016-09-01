function aggregateStimData(experiments,retinas,drugs,middleVolts)
% function [concData,timeData,voltData,infoData,ns] = aggregateStimData(experiments,retinas,drugs,stimOffsets,middleVolts)
    if nargin < 3
        stimOffsets = zeros(size(experiments));
    end
    
    nExperiments = numel(experiments);
    nDrugs = numel(drugs);
    
    if nargin < 4
        middleVolts = 8*ones(nExperiments,1);
    end
    
    uniqueDrugs = unique(drugs);
    drugExptIndices = cellfun(@(s) find(strcmp(s,drugs)),uniqueDrugs,'UniformOutput',false);
    exptsPerDrug = cellfun(@numel,drugExptIndices);
    
    ns = arrayfun(@(n) zeros(n,6),exptsPerDrug,'UniformOutput',false);
    
    nSpikesData = arrayfun(@(n) zeros(n,6,6),exptsPerDrug,'UniformOutput',false);
    dSpikesData = arrayfun(@(n) zeros(n,6,6),exptsPerDrug,'UniformOutput',false);
    snrData = arrayfun(@(n) zeros(n,6,6),exptsPerDrug,'UniformOutput',false);
    snrFinData = arrayfun(@(n) zeros(n,6,6),exptsPerDrug,'UniformOutput',false);
    threshData = arrayfun(@(n) zeros(n,6),exptsPerDrug,'UniformOutput',false);
    typesData = arrayfun(@(n) zeros(n,3,6),exptsPerDrug,'UniformOutput',false);
    
    sdData = arrayfun(@(n) zeros(n,2),exptsPerDrug,'UniformOutput',false);
    eccentricityData = arrayfun(@(n) zeros(n,2),exptsPerDrug,'UniformOutput',false);
    aspectRatioData = arrayfun(@(n) zeros(n,2),exptsPerDrug,'UniformOutput',false);
    
    infoData = arrayfun(@(n) zeros(n,2),exptsPerDrug,'UniformOutput',false);
    allInfoData = arrayfun(@(n) zeros(n,2,2),exptsPerDrug,'UniformOutput',false);
    modData = arrayfun(@(n) zeros(n,2),exptsPerDrug,'UniformOutput',false);
    perfData = arrayfun(@(n) zeros(n,2,2),exptsPerDrug,'UniformOutput',false);
    
    nUniqueDrugs = numel(uniqueDrugs);
    
    nPs = zeros(nUniqueDrugs,1);
    snrPs = zeros(nUniqueDrugs,1);
    threshPs = zeros(nUniqueDrugs,1);
    sdPs = zeros(nUniqueDrugs,1);
    eccentricityPs = zeros(nUniqueDrugs,1);
    barPs = zeros(nUniqueDrugs,1);
    
    for ii = 1:nUniqueDrugs;
        missingFF = 0;
        missingSTA = 0;
        missingBars = 0;
        
        for jj = 1:exptsPerDrug(ii)
            exptIndex = drugExptIndices{ii}(jj);
            exptDir = sprintf('JBOG%04d',experiments(exptIndex));
            
            ffFile = sprintf('%s\\Retina %c_uled_responsive_cells.mat',exptDir,retinas{exptIndex});
            
            if ~exist(ffFile,'file')
                error('Missing full-field responses file for experiment %d (%s Retina %c)',exptIndex,exptDir,retinas{exptIndex});
            end
            
            ffIndices = [3 6 7 8 11 14]-(8-middleVolts(exptIndex))*[1 0 0 0 1 0];
            
            load(ffFile);
            
            ns{ii}(jj,:) = cellfun(@numel,responsive(ffIndices)); %#ok<USENS>
            
            responseIndices = intersect(responsive{ffIndices([1 5])});
            nResponses = numel(responseIndices);
            
            if nResponses == 0
                warning('No common responders between drug and control for experiment %d (%s Retina %c)',exptIndex,exptDir,retinas{exptIndex});
                nSpikesData{ii}(jj-missingFF,:,:) = [];
                dSpikesData{ii}(jj-missingFF,:,:) = [];
                snrData{ii}(jj-missingFF,:,:) = [];
                snrFinData{ii}(jj-missingFF,:,:) = [];
                threshData{ii}(jj-missingFF,:) = [];
                missingFF = missingFF + 1;
            else
                mSpikes = median(nSpikes(:,:,responseIndices,ffIndices),1); %#ok<NODEF>
                nSpikesData{ii}(jj-missingFF,:,:) = median(mSpikes,3);
                dSpikesData{ii}(jj-missingFF,:,:) = median(mSpikes-repmat(reshape(meds(responseIndices,ffIndices),[1 1 nResponses 6]),[1 6 1 1]),3);
                snrData{ii}(jj-missingFF,:,:) = median(snr(:,responseIndices,ffIndices),2); %#ok<NODEF>
                
                for kk = 1:6
                    for ll = 1:numel(ffIndices)
                        s = snr(kk,responseIndices,ffIndices(ll));
                        snrFinData{ii}(jj-missingFF,kk,ll) = median(s(isfinite(s)));
                    end
                end
                
                threshData{ii}(jj-missingFF,:) = median(thresholds(responseIndices,ffIndices),1);
            end
            
%             allResponders = [];
%             
%             for kk = 1:6
%                 allResponders = union(allResponders,responsive{kk});
%             end
%             
%             types = zeros(3,6);
%             
%             for kk = 1:6
%                 types(1,kk) = numel(intersect(responsive{3},responsive{idx(kk)}));
%                 types(2,kk) = numel(setdiff(responsive{3},responsive{idx(kk)}));
%                 types(3,kk) = numel(setdiff(responsive{idx(kk)},responsive{3}));
%             end
%             
%             types = types./repmat(sum(types),3,1);
            
            staFile = sprintf('%s\\Retina %c_uled_flashing_squaressta.mat',exptDir,retinas{exptIndex});
            
            if ~exist(staFile,'file')
                warning('Missing STA file for experiment %d (%s Retina %c)',exptIndex,exptDir,retinas{exptIndex});
                sdData{ii}(jj-missingSTA,:) = [];
                eccentricityData{ii}(jj-missingSTA,:) = [];
                aspectRatioData{ii}(jj-missingSTA,:) = [];
                missingSTA = missingSTA + 1;
            else
                load(staFile);
                
                sdXY = cat(3,sdX,sdY);
                sd = geomean(sdXY,3); % geometric mean

                sdData{ii}(jj-missingSTA,:) = median(sd(responseIndices,:),1);
                
                a = max(sdXY,[],3);
                b = min(sdXY,[],3);
                
                eccentricity = sqrt(1-(b./a).^2);
                eccentricity(isnan(eccentricity)) = 0; % assume infinitely wide Gaussians are circular
                
                eccentricityData{ii}(jj-missingSTA,:) = median(eccentricity(responseIndices,:),1);
                
                aspectRatio = a./b;
                aspectRatio(isnan(aspectRatio)) = 1; % assume infinitely wide Gaussians are circular
                
                aspectRatioData{ii}(jj-missingSTA,:) = median(aspectRatio(responseIndices,:),1);
            end
            
            barFile = sprintf('%s\\Retina %c_bar_responses_mbg_allresp.mat',exptDir,retinas{exptIndex});
            
            if ~exist(barFile,'file')
                warning('Missing bar responses file for experiment %d (%s Retina %c)',exptIndex,exptDir,retinas{exptIndex});
                infoData{ii}(jj-missingBars,:,:) = [];
                allInfoData{ii}(jj-missingBars,:,:) = [];
                modData{ii}(jj-missingBars,:,:) = [];
                perfData{ii}(jj-missingBars,:,:) = [];
                missingBars = missingBars + 1;
                continue;
            end
            
            load(barFile);
            
%             mutualInformation = cat(3,mutualInformation{:});
            
            infoData{ii}(jj-missingBars,:) = cellfun(@(mi) median(diff(mi,[],2),1),mutualInformation);
            allInfoData{ii}(jj-missingBars,:,:) = allMutualInformation;
            modData{ii}(jj-missingBars,:) = cellfun(@(md) median(diff(md,[],2),1),modulationDepth);
            perfData{ii}(jj-missingBars,:,:) = decoderPerformance;
        end
        
        n = size(snrData{ii},1); %exptsPerDrug(ii);
        nPs(ii) = friedman(ns{ii},1,'off');
        snrPs(ii) = friedman(reshape(snrData{ii},n*6,6),n,'off');
        threshPs(ii) = friedman(threshData{ii},1,'off');
        
        sdPs(ii) = ranksum(sdData{ii}(:,1),sdData{ii}(:,2));
        eccentricityPs(ii) = ranksum(eccentricityData{ii}(:,1),eccentricityData{ii}(:,2));
        
        barPs(ii) = friedman([perfData{ii}(:,:,1); perfData{ii}(:,:,2)],size(perfData{ii},1),'off');
    end
    
    save('all_stim_data_mbg_all_resp_full_prior.mat','infoData','allInfoData','perfData','modData','ns','nSpikesData','dSpikesData','snrData','snrFinData','threshData','sdData','eccentricityData','aspectRatioData','nPs','snrPs','threshPs','sdPs','eccentricityPs','barPs');
    
    return;
    
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
        exptDir = sprintf('JBOG%04d',experiments(ii));
        responsiveCellss = cell(nDrugs,1);
        
        for jj = 1:nDrugs
            responsiveCellFile = sprintf('%s\\%s_%s_responsive_cells.mat',exptDir,exptDir,drugs{jj});
            
            if exist(responsiveCellFile,'file')
                responsiveCells = load(responsiveCellFile,'responsiveCells');
                responsiveCellss{jj} = responsiveCells.responsiveCells;
                nCells = size(responsiveCells,1);
                data{ii,jj} = zeros(nCells,6,5,6,4);
            end
        end
        
        recDirs = dir(sprintf('%s\\Stim*',exptDir));
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
            
            dataFile = sprintf('%s\\%s\\%s_uled_square_responses_newmethod.mat',exptDir,recDir,recDir);
            
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
                responsiveCellFile = sprintf('%s\\%s_%s_%d_responsive_cells.mat',exptDir,exptDir,drugs{drugIndex},drugSuffix);
                
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