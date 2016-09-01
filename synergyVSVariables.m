isEvenBinning = false;
filePrefix = 'P39';

%%

clear idx idxBOTH;
load('ALLunits.mat','idx','idxBOTH')

if exist('idxBOTH','var') && ~exist('idx','var')
    idx = idxBOTH;
end

load('channelNames.mat')

if exist('./IDXalternative.mat','file')
    load('./IDXalternative.mat');
    idx = idxGood;
    figSuffix = '_alt';
end

% SelectN = dir('*SelectN.mat');
% load(SelectN.name);
% 
% load('synergy2.mat')

if exist('./synergy25.mat','file')
    load('synergy25.mat')
else
    load([filePrefix '_synergy25.mat'])
end

%%

% n = numel(SelectN);
n = numel(idx);
N = n/2*(n-1);

valid = [channelNames{6,:}] == 0;

assert(all(ismember(idx-1,find(valid))),'Selected units should only included non-marked units');

indices = zeros(n,1);

for ii = 1:n
%     indices(ii) = find(ismember(find(valid),SelectN(ii)));
    indices(ii) = find(ismember(find(valid),idx(ii)-1));
end

%%

pairs = zeros(N,2);

nn = 0;
for ii = 1:n-1
    for jj = (ii+1):n
        nn = nn + 1;
        pairs(nn,:) = [ii jj];
    end
end

%%

% OFF = biasIndex <= -0.25;
% 
% t2p = [ALLunits{SelectN+1,16}]';
% fast = t2p < 150;
% 
% fastOFF = find(fast & OFF);
% 
% fastOFFPairs = ismember(pairs(:,1),fastOFF) & ismember(pairs(:,2),fastOFF);
% pairs = pairs(fastOFFPairs,:);

%%

% figure;
% barwitherr(squeeze(std(pids(fastOFFPairs,:,:)))/sqrt(sum(fastOFFPairs)),squeeze(mean(pids(fastOFFPairs,:,:))));
% set(gca,'XTickLabel',{'Redundancy' 'Unique 1' 'Unique 2' 'Synergy' 'Information'})
% ylabel(ylabels{3});
% close(gcf);

%%

% y = double(vertcat(channelNames{2,SelectN}));
% x = double(vertcat(channelNames{3,SelectN}));
y = double(vertcat(channelNames{2,idx-1}));
x = double(vertcat(channelNames{3,idx-1}));

dx = x(pairs(:,2))-x(pairs(:,1));
dy = y(pairs(:,2))-y(pairs(:,1));

dd = sqrt(dx.^2+dy.^2);
dt = atan2(dy,dx);

biasIndex = vertcat(ALLunits{idx,11});
% biasIndex = vertcat(ALLunits{SelectN+1,12});
db = abs(biasIndex(pairs(:,2))-biasIndex(pairs(:,1)));

%%

titles = {'Redundancy' 'Synergy' 'Total Information' 'Number of Pairs'};
ylabels = {'Proportion of Information' 'Proportion of Information' 'Information (bits)' '# Pairs'};

legends = cell(4,1);
sfs = [75 37.5 18.4 9.4];
for ii = 1:4
    legends{ii} = sprintf('%3.1f mcpd',sfs(ii));
end

%%

dns = zeros(89,1);
diqrs = zeros(89,3,4,3);

if isEvenBinning
    dedges = 0:89; %#ok<*UNRCH>
else
    dedges = prctile(dd,linspace(0,100,90));
end

for ii = 1:89
    if ii == 89
        dataIndices = dd >= dedges(ii) & dd <= dedges(ii+1);
    else
        dataIndices = dd >= dedges(ii) & dd < dedges(ii+1);
    end
    
    dns(ii) = sum(dataIndices);
    
    info = pids(dataIndices,5,:);
    diqrs(ii,3,:,:) = squeeze(prctile(info,[25 50 75]))';
    
    data = pids(dataIndices,[1 4],:);
    data = data./repmat(info,[1 2 1]);
    data(isnan(data)) = 0;
    diqrs(ii,1:2,:,:) = permute(prctile(data,[25 50 75]),[2 3 1]);
end

%%

dcentres = mean([dedges(1:end-1); dedges(2:end)])*0.042;

figure;
set(gcf,'Position',[1 49 1920 946]);
    
for ii = 1:3
    subplot(2,2,ii)
    boundedline(dcentres,squeeze(diqrs(:,ii,:,2)),permute(squeeze(diqrs(:,ii,:,[2 3])-diqrs(:,ii,:,[1 2])),[1 3 2]),'alpha');
end

subplot(2,2,4);
plot(dcentres,dns);

for ii = 1:4
    subplot(2,2,ii)
    title(titles{ii});
    ylabel(ylabels{ii});
    xlim([0 2.76]);
end

subplot(2,2,3);
legend(legends,'Location','NorthWest')
xlabel('Distance (mm)')

figFile = sprintf('%s_synergy25vsdist_long',filePrefix);
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');

%%

for ii = 1:4
    subplot(2,2,ii);
    xlim([0 1]);
end

figFile = sprintf('%s_synergy25vsdist_short',filePrefix);
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');

close(gcf);

%%

tns = zeros(100,1);
tiqrs = zeros(100,3,4,3);

if isEvenBinning
    tedges = linspace(0,pi,101);
else
    tedges = prctile(dt,0:100);
end

badBins = find(diff(tedges) == 0);
badStarts = badBins([1 find(diff(badBins) > 1)+1]);
badEnds = badBins([find(diff(badBins) > 1) end])+1;

%%

if ~isempty(badBins)
    for ii = 1:numel(badStarts)
        if badStarts(ii) == 1
            startIndex = 1;
        else
            startIndex = badStarts(ii)-1;
        end
        
        if badEnds(ii) == numel(tedges)
            endIndex = numel(tedges);
        else
            endIndex = badEnds(ii)+1;
        end
        
        tedges(startIndex:endIndex) = linspace(tedges(startIndex),tedges(endIndex),endIndex-startIndex+1);
    end
end

%%

for ii = 1:100
    if ii == 100
        dataIndices = dt >= tedges(ii) & dt <= tedges(ii+1);
    else
        dataIndices = dt >= tedges(ii) & dt < tedges(ii+1);
    end
    
    tns(ii) = sum(dataIndices);
    
    info = pids(dataIndices,5,:);
    tiqrs(ii,3,:,:) = squeeze(prctile(info,[25 50 75]))';
    
    data = pids(dataIndices,[1 4],:);
    data = data./repmat(info,[1 2 1]);
    data(isnan(data)) = 0;
    tiqrs(ii,1:2,:,:) = permute(prctile(data,[25 50 75]),[2 3 1]);
end

%%

figure;
set(gcf,'Position',[1 49 1920 946]);

tcentres = mean([tedges(1:end-1); tedges(2:end)]);

for ii = 1:3
    subplot(2,2,ii)
    boundedline(tcentres,squeeze(tiqrs(:,ii,:,2)),permute(squeeze(tiqrs(:,ii,:,[2 3])-tiqrs(:,ii,:,[1 2])),[1 3 2]),'alpha');
end

subplot(2,2,4);
plot(tcentres,tns);

for ii = 1:4
    subplot(2,2,ii)
    title(titles{ii});
    ylabel(ylabels{ii});
    xlim([0 pi])
    set(gca,'XTick',(0:4)*pi/4,'XTickLabel',0:45:180);
end

subplot(2,2,3);
legend(legends,'Location','NorthWest')
xlabel('Angle (degrees)')

figFile = sprintf('%s_synergy25vsangle',filePrefix);
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');

close(gcf);

%%

bns = zeros(20,1);
biqrs = zeros(20,3,4,3);

if isEvenBinning
    bedges = 0:0.1:2;
else
    bedges = prctile(db,0:5:100);
end

badBins = find(diff(bedges) == 0);
badStarts = badBins([1 find(diff(badBins) > 1)+1]);
badEnds = badBins([find(diff(badBins) > 1) end])+1;

%%

if ~isempty(badBins)
    for ii = 1:numel(badStarts)
        if badStarts(ii) == 1
            startIndex = 1;
        else
            startIndex = badStarts(ii)-1;
        end
        
        if badEnds(ii) == numel(bedges)
            endIndex = numel(bedges);
        else
            endIndex = badEnds(ii)+1;
        end
        
        bedges(startIndex:endIndex) = linspace(bedges(startIndex),bedges(endIndex),endIndex-startIndex+1);
    end
end

%%

for ii = 1:20
    if ii < 20
        dataIndices = db >= bedges(ii) & db < bedges(ii+1);
    else
        dataIndices = db >= bedges(ii) & db <= bedges(ii+1);
    end
    
    bns(ii) = sum(dataIndices);
    
    info = pids(dataIndices,5,:);
    biqrs(ii,3,:,:) = squeeze(prctile(info,[25 50 75]))';
    
    data = pids(dataIndices,[1 4],:);
    data = data./repmat(info,[1 2 1]);
    data(isnan(data)) = 0;
    biqrs(ii,1:2,:,:) = permute(prctile(data,[25 50 75]),[2 3 1]);
end

%%

bcentres = mean([bedges(1:end-1); bedges(2:end)]);

figure;
set(gcf,'Position',[1 49 1920 946]);

for ii = 1:3
    subplot(2,2,ii)
    boundedline(bcentres,squeeze(biqrs(:,ii,:,2)),permute(squeeze(biqrs(:,ii,:,[2 3])-biqrs(:,ii,:,[1 2])),[1 3 2]),'alpha');
end

subplot(2,2,4);
plot(bcentres,bns);

for ii = 1:4
    subplot(2,2,ii)
    title(titles{ii});
    ylabel(ylabels{ii});
    xlim([0 2]);
end

subplot(2,2,3);
legend(legends,'Location','NorthWest')
xlabel('{\Delta}Bias Index')

figFile = sprintf('%s_synergy25vsbias',filePrefix);
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');

close(gcf);

%%

dbiqrs = zeros(20,20,3,4,3);
dbns = zeros(20,20,1);

dedges = prctile(dd,linspace(0,100,21));
bedges = prctile(db,linspace(0,100,21));

badBins = find(diff(bedges) == 0);
badStarts = badBins([1 find(diff(badBins) > 1)+1]);
badEnds = badBins([find(diff(badBins) > 1) end])+1;

%%

if ~isempty(badBins)
    for ii = 1:numel(badStarts)
        if badStarts(ii) == 1
            startIndex = 1;
        else
            startIndex = badStarts(ii)-1;
        end
        
        if badEnds(ii) == numel(bedges)
            endIndex = numel(bedges);
        else
            endIndex = badEnds(ii)+1;
        end
        
        bedges(startIndex:endIndex) = linspace(bedges(startIndex),bedges(endIndex),endIndex-startIndex+1);
    end
end

%%

for ii = 1:20
    for jj = 1:20
        if ii < 20
            dIndices = dd >= dedges(ii) & dd < dedges(ii+1);
        else
            dIndices = dd >= dedges(ii) & dd <= dedges(ii+1);
        end
        
        if jj < 20
            bIndices = db >= bedges(jj) & db < bedges(jj+1);
        else
            bIndices = db >= bedges(jj) & db <= bedges(jj+1);
        end

        dataIndices = dIndices & bIndices;

        dbns(ii,jj) = sum(dataIndices);

        info = pids(dataIndices,5,:);
        dbiqrs(ii,jj,3,:,:) = squeeze(prctile(info,[25 50 75],1))';

        data = pids(dataIndices,[1 4],:);
        data = data./repmat(info,[1 2 1]);
        data(isnan(data)) = 0;
        dbiqrs(ii,jj,1:2,:,:) = permute(prctile(data,[25 50 75],1),[2 3 1]);
    end
end

%%

for ii = 1:4
    figure;
    set(gcf,'Position',[1 49 1920 946]);
    
    Cs = cell(4,1);
    for jj = 1:4
        if jj == 4
            C = dbns;
        else
            C = dbiqrs(:,:,jj,ii,2);
        end
        
        Cs{jj} = C;

        minC = min(C(:));
        maxC = max(C(:));
        subplot(2,2,jj);
        pcolor(bedges,dedges*0.042,[C zeros(20,1); zeros(1,21)]);
        caxis([minC maxC]);
        xlim(bedges([1 end]));
        ylim(dedges([1 end])*0.042);
        xlabel('{\Delta}Bias Index');
        ylabel('Distance (mm)');
        colorbar;
        title(titles{jj});
        shading flat;
    end
    
    figPrefix = strrep(sprintf('%s_synergy25_vs_bias_dist_%3.1fmcpd',filePrefix,sfs(ii)),'.','_');
    figFile = [figPrefix '_allbins'];
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');

    close(gcf);
    
    continue;
    
    for jj = 1:4
        C = Cs{jj}(1:end-2,1:end-4);
        minC = min(C(:));
        maxC = max(C(:));
        subplot(2,2,jj);
        caxis([minC maxC]);
        xlim(bedges([1 end-4]));
        ylim(dedges([1 end-2])*0.042);
        colorbar;
    end
    
    figFile = [figPrefix '_notlastbins'];
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
%     suptitle(legends{ii});
end