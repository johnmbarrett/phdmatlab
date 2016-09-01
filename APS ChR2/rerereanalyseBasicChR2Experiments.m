% P = 95;
P = 96;
% P = 49;
% P = 52;

if isnumeric(P)
    prefix = sprintf('P%d',P);
else
    prefix = P;
end

%%

tic;
if exist('detectionFile','var')
    if ischar(detectionFile)
        load(detectionFile,'allBestUnits');
    else
        uniqueDetectionFile = unique(detectionFile);
        
        if numel(uniqueDetectionFile) == 1
            load(uniqueDetectionFile{1},'allBestUnits');
        else
            allAllBestUnits = cell(numel(detectionFile),2);

            for ii = 1:numel(detectionFile)
                load(detectionFile{ii},'allBestUnits');

                allAllBestUnits(ii,:) = allBestUnits;
            end

            allBestUnits = allAllBestUnits;
        end
    end
else
    load('ChR2_cells_new');
end

load('EventNo');
load('spiketimestamps.mat');
toc;

tic;
if isnumeric(P) && P == 96
    recordings = [1 2 3 5];
elseif ~exist('recordings','var')
    recordings = 1:4;
end

% EventNo = EventNo(recordings);

% fin = fopen('originalSplit.txt');
% closeFile = onCleanup(@() fclose(fin));
% nSamples = textscan(fin, '%*s %d', 'CommentStyle', '#', 'Delimiter', ',');
% sampleRate = 7055.258405732180;
% endTimes = cumsum(double(nSamples{1}))/sampleRate;
toc;

%%

tic;
% n = numel(bestUnits);
nRecs = size(allBestUnits,1);
afterFR = cell(nRecs,2);
% k = 5;
beforeFR = cell(nRecs,2);

if ~exist('stimIndices','var')
    stimIndices = repmat({1:2:59 2:2:60},nRecs,1);
end

for ii = 1:nRecs
    for jj = 1:2
        afterFR{ii,jj} = cellfun(@(s) mean(arrayfun(@(t) sum(s > t & s <= t+0.25),EventNo{recordings(ii)}(stimIndices{ii,jj}(2:2:end))))/0.25,spiketimestamps(allBestUnits{ii,jj}))';
        beforeFR{ii,jj} = cellfun(@(s) mean(arrayfun(@(t) sum(s > t-0.5 & s <= t),EventNo{recordings(ii)}(stimIndices{ii,jj}(2:2:end))))/0.5,spiketimestamps(allBestUnits{ii,jj}))';
%         afterFR{ii,jj} = cellfun(@(S) (sum(mean(S{jj}(2:2:end,1:2),1),2)+mean(S{jj}(2:2:end,3))/2)/0.25,nSpikesss(ismember(goodUnits,allBestUnits{ii,jj}),ii));
%         beforeFR{ii,jj} = cellfun(@(S) sum(mean(S{3-jj}(2:2:end,16:end),1),2)/0.5,nSpikesss(ismember(goodUnits,allBestUnits{ii,jj}),ii));
    end
end
toc;

% bad = cellfun(@(p) p < 0.05,pSpont,'UniformOutput',false); %min(spontFR,[],2) == 0;
% disp(n);
% disp(sum(bad));
% disp(n-sum(bad));

%%

% if P == 96
%     pspont = prctile(spontFR,[25 50 75]);
%     mspont = squeeze(pspont(2,:))';
%     uspont = squeeze(pspont(3,:))'-mspont;
%     lspont = mspont-squeeze(pspont(1,:))';
%     
%     figure;
%     barwitherr([lspont uspont],mspont);
% end

%%

tic;
% evokedPFR = cellfun(@(E,S) E./S,evokedFR,spontFR,'UniformOutput',false);
evokedDFR = cellfun(@(A,B) A-B,afterFR,beforeFR,'UniformOutput',false);
evokedADFR = cellfun(@(A,B) abs(A-B),afterFR,beforeFR,'UniformOutput',false);
evokedPFR = cellfun(@(A,B) 100*(A./B-1),afterFR,beforeFR,'UniformOutput',false);
evokedAPFR = cellfun(@(A,B) 100*abs(A./B-1),afterFR,beforeFR,'UniformOutput',false);
% evokedPFR = cellfun(@(E,S,bad) E(~bad)./S(~bad),evokedFR,repmat(spontFR,1,2),repmat(bad,1,2),'UniformOutput',false);
% peakDFR = peakDFR(~bad,:,:);

%%

pdfr = cell2mat(reshape(cellfun(@(E) prctile(E,[25;50;75]),evokedDFR,'UniformOutput',false),[1 nRecs 2]));
mdfr = squeeze(pdfr(2,:,:));
udfr = squeeze(pdfr(3,:,:))-mdfr;
ldfr = mdfr-squeeze(pdfr(1,:,:));
toc;

%%

tic;
if isnumeric(P)
    if P == 49
        flab = 'Blockers';
        xlab = 'Neutral Density Filter';
        xticks = [4.5 4.5 2.2 2.2];
        xfill = [1.5 3.5; 1.5 3.5; 2.5 5; 2.5 5];
        yfill = [0.1 40 40 0.1; 0.1 40 40 0.1]';
        ylims = [0 25];
    elseif P == 52 || P == 95 
        flab = 'Blockers';
        xlab = 'Neutral Density Filter';
        xticks = [4.5 2.2 4.5 2.2];
        xfill = [2.5 2.5 5 5];
        yfill = [0.1 40 40 0.1];
        ylims = [0 25];
    elseif P == 96
        flab = 'MFA';
        xlab = 'KCl Concentration (mM)';
        xticks = [3 3 6 9];
        xfill = [1.5 1.5 5 5];
        yfill = [0.1 300 300 0.1];
        ylims = [0 200];
    end
end
toc;

%%

tic;
figure;

hold on;
hf = fill(xfill,yfill,[0.75 1 0.75],'EdgeColor','none');
hb = barwitherr(cat(3,ldfr,udfr),mdfr);
set(get(gca,'Children'),'LineWidth',1.5);
set(gca,'LineWidth',1.5,'XTick',1:4,'XTickLabel',xticks);
xlabel(xlab);
xlim([0.5 4.5]);
ylabel('Change in Firing Rate');
ylim([1.1*min(pdfr(:)) 1.1*max(pdfr(:))]);

legend([hb hf(1)],{'ON Responses' 'OFF Responses' flab},'LineWidth',1.5,'Location','NorthWest');
toc;

%%

% figure;
% distributionPlot(reshape(evokedDFR',8,1),'color',repmat({'b';'r'},4,1),'histOpt',2,'showMM',6)

%%

tic;
figName = sprintf('%s_responses_new2_250ms',prefix);
saveas(gcf,figName,'fig');
export_fig(figName,'-eps','-png','-transparent','-painters');
toc;

%%

save(sprintf('%s_responses_new.mat',prefix),'-v7.3','afterFR','beforeFR','evokedDFR','evokedADFR','evokedPFR','evokedAPFR','pdfr','mdfr','ldfr','udfr');