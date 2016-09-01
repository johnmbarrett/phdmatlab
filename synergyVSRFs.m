isEvenBinning = false;
filePrefix = 'P39';

%%

load('channelNames.mat')

if exist([filePrefix '_synergy25_alt.mat'],'file')
    load([filePrefix '_synergy25_alt.mat']);
elseif exist([filePrefix '_synergy25.mat'],'file')
    load([filePrefix '_synergy25.mat']);
elseif exist('synergy25.mat','file')
    load('synergy25.mat');
else
    error('Couldn''t find synergy file');
end

%%

if exist('IDXalternative.mat','file')
    load('ALLunits25.mat','ALLunits');
    load('IDXalternative.mat');
    auidx = idxGood-1;
else
    clear idx idxBOTH
    load('ALLunits25.mat','ALLunits','idxBOTH','idx')

    if exist('idxBOTH','var') && ~exist('idx','var')
        auidx = idxBOTH-1;
    else
        auidx = idx-1;
    end
end

%%

rffile = dir('RF_selected*.mat');

if numel(rffile) ~= 1
    return;
end

load(rffile.name);

%%

RFs = vertcat(RF_selected{:});
rfidx = vertcat(RFs.id);
RFs = RFs(ismember(rfidx,auidx));
rfidx = rfidx(ismember(rfidx,auidx));

%%

n = numel(rfidx);
N = n/2*(n-1);

valid = find([channelNames{6,:}] == 0);

m = numel(auidx);
M = m/2*(m-1);

assert(all(ismember(rfidx,valid)),'Selected units should only included non-marked units');
assert(M == size(pids,1),'Number of selected unit pairs doesn''t match number of rows in PID, probably loaded the wrong synergy file');

%%

pairs = zeros(N,2);

nn = 0;
for ii = 1:n-1
    for jj = (ii+1):n
        nn = nn + 1;
        pairs(nn,:) = [ii jj];
    end
end

pairIndices = zeros(N,1);

nn = 0;
mm = 0;
for ii = 1:m-1
    for jj = (ii+1):m
        mm = mm + 1;
        
        if ismember(auidx(ii),rfidx) && ismember(auidx(jj),rfidx)
            nn = nn + 1;
            pairIndices(nn) = mm;
        end
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

xy = vertcat(RFs.Centroid);
dxy = (xy(pairs(:,1),:)-xy(pairs(:,2),:))*(664/27)*4;

dd = sqrt(sum(dxy.^2,2));
dt = atan(dxy(:,2)./dxy(:,1)); % NOT atan2, because the angle shouldn't depend on the pair ordering

biasIndex = vertcat(ALLunits{rfidx+1,11});
db = abs(biasIndex(pairs(:,2))-biasIndex(pairs(:,1)));

%%

titles = {'Redundancy' 'Synergy' 'Total Information' 'Number of Pairs'};
ylabels = {'Proportion of Information' 'Proportion of Information' 'Information (bits)' '# Pairs'};

legends = cell(4,1);
% sfs = [75 37.5 18.4 9.4];
bws = 100*2.^(1:4);
for ii = 1:4
%     legends{ii} = sprintf('%3.1f mcpd',sfs(ii));
    legends{ii} = sprintf('%d {\\mu}m',bws(ii));
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
    
    info = pids(pairIndices(dataIndices),5,:);
    diqrs(ii,3,:,:) = squeeze(prctile(info,[25 50 75]))';
    
    data = pids(pairIndices(dataIndices),[1 4],:);
    data = data./repmat(info,[1 2 1]);
    data(isnan(data)) = 0;
    diqrs(ii,1:2,:,:) = permute(prctile(data,[25 50 75]),[2 3 1]);
end

%%

dcentres = mean([dedges(1:end-1); dedges(2:end)]); %*0.042;

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
    xlim(dcentres([1 end]));
end

subplot(2,2,3);
legend(legends,'Location','NorthWest')
xlabel('Distance ({\mu}m)')

figFile = sprintf('%s_synergy25vsrfdist_long',filePrefix);
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');
% close(gcf);

%%

for ii = 1:4
    subplot(2,2,ii);
    xlim([0 500]);
end

figFile = sprintf('%s_synergy25vsrfdist_short',filePrefix);
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
tedges(badBins+1) = mean([tedges(badBins); tedges(badBins+2)]);

for ii = 1:100
    if ii == 100
        dataIndices = dt >= tedges(ii) & dt <= tedges(ii+1);
    else
        dataIndices = dt >= tedges(ii) & dt < tedges(ii+1);
    end
    
    tns(ii) = sum(dataIndices);
    
    info = pids(pairIndices(dataIndices),5,:);
    tiqrs(ii,3,:,:) = squeeze(prctile(info,[25 50 75]))';
    
    data = pids(pairIndices(dataIndices),[1 4],:);
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
    xlim([-pi/2 pi/2])
    set(gca,'XTick',(-2:2)*pi/4,'XTickLabel',-90:45:90);
end

subplot(2,2,3);
legend(legends,'Location','NorthWest')
xlabel('Angle Relative to Gratings (degrees)')

figFile = sprintf('%s_synergy25vsrfangle',filePrefix);
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
    
    info = pids(pairIndices(dataIndices),5,:);
    biqrs(ii,3,:,:) = squeeze(prctile(info,[25 50 75]))';
    
    data = pids(pairIndices(dataIndices),[1 4],:);
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

figFile = sprintf('%s_synergy25vsrfbias',filePrefix);
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');

close(gcf);

%%

centres = vertcat(RFs.Centroid);
x0 = centres(:,1);
y0 = centres(:,2);
A = vertcat(RFs.MajorAxisLength);
B = vertcat(RFs.MinorAxisLength);
phi = vertcat(RFs.Orientation);

%%

alphax = atan2(B.*sin(phi),A.*cos(phi));
alphay = atan2(B.*cos(phi),A.*sin(phi));

thetax = [pi-alphax -alphax];
thetay = [pi+alphay  alphay];

xlims = [x0 x0]+[A A].*cos(thetax).*cos([phi phi])-[B B].*sin(thetax).*sin([phi phi]);
ylims = [y0 y0]+[A A].*cos(thetay).*sin([phi phi])+[B B].*sin(thetay).*cos([phi phi]);

% TODO : is this true?
assert(all(xlims(:,1) < xlims(:,2)) && all(ylims(:,1) < ylims(:,2)));

%%

isPossibleOverlap = find( ...
      (xlims(pairs(:,1),1) < xlims(pairs(:,2),2)) & (xlims(pairs(:,1),2) > xlims(pairs(:,2),1)) ...
    & (ylims(pairs(:,1),1) < ylims(pairs(:,2),2)) & (ylims(pairs(:,1),2) > ylims(pairs(:,2),1)) ...
    );

o = zeros(N,1);

%%

No = numel(isPossibleOverlap);

xmin = min(xlims(pairs(isPossibleOverlap,1),1),xlims(pairs(isPossibleOverlap,2),1));
xmax = max(xlims(pairs(isPossibleOverlap,1),2),xlims(pairs(isPossibleOverlap,2),2));
ymin = min(ylims(pairs(isPossibleOverlap,1),1),ylims(pairs(isPossibleOverlap,2),1));
ymax = max(ylims(pairs(isPossibleOverlap,1),2),ylims(pairs(isPossibleOverlap,2),2));

Nxy = 1e5;
rx = xmax-xmin;
ry = ymax-ymin;
delta = (-(rx+ry)+sqrt((rx+ry).^2+4*rx.*ry*Nxy))./(2*rx.*ry);

%%

for ii = 1:No
    idx1 = pairs(isPossibleOverlap(ii),1);
    idx2 = pairs(isPossibleOverlap(ii),2);
    
    tic;
    [Y,X] = ndgrid(ymin(ii):1/delta(ii):ymax(ii),xmin(ii):1/delta(ii):xmax(ii));
    n(ii,jj) = numel(Y);

    E1 = ((X-x0(idx1))*cos(phi(idx1))+(Y-y0(idx1))*sin(phi(idx1))).^2/A(idx1)^2+((Y-y0(idx1))*cos(phi(idx1))-(X-x0(idx1))*sin(phi(idx1))).^2/B(idx1)^2 <= 1;
    E2 = ((X-x0(idx2))*cos(phi(idx2))+(Y-y0(idx2))*sin(phi(idx2))).^2/A(idx2)^2+((Y-y0(idx2))*cos(phi(idx2))-(X-x0(idx2))*sin(phi(idx2))).^2/B(idx2)^2 <= 1;
    o(isPossibleOverlap(ii)) = sum(sum(E1 & E2))/sum(sum(E1 | E2));
    fprintf('Pair %d, points %d, time %f seconds\n',ii,numel(Y),toc);
end

%%

isOverlapping = find(o > 0);
meanInfoVOverlap = [mean(pids(pairIndices(isOverlapping),:,:)); mean(pids(pairIndices(o == 0),:,:))];
stdInfoVOverlap = [std(pids(pairIndices(isOverlapping),:,:)); std(pids(pairIndices(o == 0),:,:))];

%%

figure;
set(gcf,'Position',[1 49 1920 946]);

for ii = 1:4
    subplot(2,2,ii);
    barwitherr(stdInfoVOverlap(:,:,ii)'/sqrt(sum(isOverlapping)),meanInfoVOverlap(:,:,ii)');
    set(gca,'XTickLabel',{'Redundancy' 'Unique 1' 'Unique 2' 'Synergy' 'Total Information'});
    title(sprintf('%3.1f mcpd',sfs(ii)));
end

legend({'Overlapping' 'Non-overlapping'},'Location','NorthWest');

subplot(2,2,3);
ylabel('Information');

figFile = sprintf('%s_synergy25vsisoverlapping',filePrefix);
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');

close(gcf);

%%

overlap = o(o > 0);
ons = zeros(20,1);
oiqrs = zeros(20,3,4,3);

if isEvenBinning
    oedges = 0:0.05:1;
else
    oedges = prctile(overlap,0:5:100);
end

badBins = find(diff(oedges) == 0);

%%

if ~isempty(badBins)
    badStarts = badBins([1 find(diff(badBins) > 1)+1]);
    badEnds = badBins([find(diff(badBins) > 1) end])+1;
    
    for ii = 1:numel(badStarts)
        if badStarts(ii) == 1
            startIndex = 1;
        else
            startIndex = badStarts(ii)-1;
        end
        
        if badEnds(ii) == numel(oedges)
            endIndex = numel(oedges);
        else
            endIndex = badEnds(ii)+1;
        end
        
        oedges(startIndex:endIndex) = linspace(oedges(startIndex),oedges(endIndex),endIndex-startIndex+1);
    end
end

%%

for ii = 1:20
    if ii < 20
        dataIndices = overlap >= oedges(ii) & overlap < oedges(ii+1);
    else
        dataIndices = overlap >= oedges(ii) & overlap <= oedges(ii+1);
    end
    
    ons(ii) = sum(dataIndices);
    
    info = pids(pairIndices(isOverlapping(dataIndices)),5,:);
    oiqrs(ii,3,:,:) = squeeze(prctile(info,[25 50 75]))';
    
    data = pids(pairIndices(isOverlapping(dataIndices)),[1 4],:);
    data = data./repmat(info,[1 2 1]);
    data(isnan(data)) = 0;
    oiqrs(ii,1:2,:,:) = permute(prctile(data,[25 50 75]),[2 3 1]);
end

%%

ocentres = mean([oedges(1:end-1); oedges(2:end)]);

figure;
set(gcf,'Position',[1 49 1920 946]);

for ii = 1:3
    subplot(2,2,ii)
    boundedline(ocentres,squeeze(oiqrs(:,ii,:,2)),permute(squeeze(oiqrs(:,ii,:,[2 3])-oiqrs(:,ii,:,[1 2])),[1 3 2]),'alpha');
end

subplot(2,2,4);
plot(ocentres,ons);

for ii = 1:4
    subplot(2,2,ii)
    title(titles{ii});
    ylabel(ylabels{ii});
    xlim([0 1]);
end

subplot(2,2,3);
legend(legends,'Location','NorthWest')
xlabel('RF Overlap')

figFile = sprintf('%s_synergy25vsrfoverlap',filePrefix);
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');

% close(gcf);

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

        info = pids(pairIndices(dataIndices),5,:);
        dbiqrs(ii,jj,3,:,:) = squeeze(prctile(info,[25 50 75],1))';

        data = pids(pairIndices(dataIndices),[1 4],:);
        data = data./repmat(info,[1 2 1]);
        data(isnan(data)) = 0;
        dbiqrs(ii,jj,1:2,:,:) = permute(prctile(data,[25 50 75],1),[2 3 1]);
    end
end

%%

figs = zeros(1,4);
    
minCs = Inf(1,4);
maxCs = -Inf(1,4);
    
for ii = 1:4
    figs(ii) = figure;
    set(gcf,'Position',[1 49 1920 946]);
    
    for jj = 1:4
        if jj == 4
            C = dbns;
        else
            C = dbiqrs(:,:,jj,ii,2);
        end
        
        minC = min(C(:));
        maxC = max(C(:));
        minCs(jj) = min(minCs(jj),minC);
        maxCs(jj) = max(maxCs(jj),maxC);
        subplot(2,2,jj);
        pcolor(bedges,dedges,[C zeros(20,1); zeros(1,21)]);
        caxis([minC maxC]);
        xlim(bedges([1 end]));
        ylim(dedges([1 end]));
        xlabel('{\Delta}Bias Index');
        ylabel('Distance ({\mu}m)');
        colorbar;
        title(titles{jj});
        shading flat;
    end
    
    figPrefix = strrep(sprintf('%s_synergy25_vs_bias_rfdist_%dum',filePrefix,bws(ii)),'.','_');
    figFile = [figPrefix '_diffcaxis'];
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
end

for ii = 1:4
    figure(figs(ii));
    
    for jj = 1:4
        subplot(2,2,jj);
        caxis([minCs(jj) maxCs(jj)]);
    end
    
%     figPrefix = strrep(sprintf('%s_synergy25_vs_bias_rfdist_%3.1fmcpd',filePrefix,sfs(ii)),'.','_');
    figPrefix = strrep(sprintf('%s_synergy25_vs_bias_rfdist_%dum',filePrefix,bws(ii)),'.','_');
    figFile = [figPrefix '_samecaxis'];
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');

    close(gcf);
end
    
%     continue;
%     
%     for jj = 1:4
%         C = Cs{jj}(1:end-2,1:end-4);
%         minC = min(C(:));
%         maxC = max(C(:));
%         subplot(2,2,jj);
%         caxis([minC maxC]);
%         xlim(bedges([1 end-4]));
%         ylim(dedges([1 end-2])*0.042);
%         colorbar;
%     end
%     
%     figFile = [figPrefix '_notlastbins'];
%     saveas(gcf,figFile,'fig');
%     saveas(gcf,figFile,'png');
% %     suptitle(legends{ii});
% end