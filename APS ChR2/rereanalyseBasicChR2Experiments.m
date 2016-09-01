% P = 95;
P = 96;

%%

load('ChR2_cells_ff');
load('EventNo');
load('spiketimestamps.mat');

if P == 96
    recordings = [1 2 3 5];
else
    recordings = 1:4;
end

EventNo = EventNo(recordings);

fin = fopen('originalSplit.txt');
closeFile = onCleanup(@() fclose(fin));
nSamples = textscan(fin, '%*s %d', 'CommentStyle', '#', 'Delimiter', ',');
sampleRate = 7055.258405732180;
endTimes = cumsum(double(nSamples{1}))/sampleRate;

%%

n = numel(bestUnits);
evokedFR = cell(4,2);
% k = 5;
spontFR = cell(4,2);
spontFREnd = cell(4,2);
pSpont = cell(4,2);

for ii = 1:4
    if ii == 1
        t0 = 0;
    else
        t0 = endTimes(recordings(ii)-1);
    end
    
    t1 = EventNo{ii}(1);
    t2 = EventNo{ii}(end);
    t3 = endTimes(recordings(ii));
    
    T = t1-t0;
    U = t3-t2;
    
    disp(T);
    disp(U);
    
    for jj = 1:2
        spontFR{ii,jj} = cellfun(@(t) sum(t > t0 & t < t1),spiketimestamps(allBestUnits{ii,jj}))'/T;
        spontFREnd{ii,jj} = cellfun(@(t) sum(t > t2 & t < t3),spiketimestamps(allBestUnits{ii,jj}))'/U;
        pSpont{ii,jj} = 1-poisscdf(spontFREnd{ii,jj},spontFR{ii,jj});
        evokedFR{ii,jj} = cellfun(@(S) sum(mean(S{jj}(2:2:end,1:k),1),2)/(k*0.1),nSpikesss(ismember(goodUnits,allBestUnits{ii,jj}),ii))-spontFR{ii,jj};
    end
end

bad = cellfun(@(p) p < 0.05,pSpont,'UniformOutput',false); %min(spontFR,[],2) == 0;
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

% evokedPFR = cellfun(@(E,S) E./S,evokedFR,spontFR,'UniformOutput',false);
% evokedPFR = cellfun(@(E,S,bad) E(~bad)./S(~bad),evokedFR,spontFR,bad,'UniformOutput',false);
% evokedPFR = cellfun(@(E,S,bad) E(~bad)./S(~bad),evokedFR,repmat(spontFR,1,2),repmat(bad,1,2),'UniformOutput',false);
% peakDFR = peakDFR(~bad,:,:);

%%

pdfr = cell2mat(reshape(cellfun(@(E,bad) prctile(E(~bad),[25;50;75]),evokedFR,bad,'UniformOutput',false),[1 4 2]));
mdfr = squeeze(pdfr(2,:,:));
udfr = squeeze(pdfr(3,:,:))-mdfr;
ldfr = mdfr-squeeze(pdfr(1,:,:));

%%

if P == 95
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

%%

figure;

hold on;
hf = fill(xfill,yfill,[0.75 1 0.75],'EdgeColor','none');
hb = barwitherr(cat(3,ldfr,udfr),mdfr);
set(get(gca,'Children'),'LineWidth',1.5);
set(gca,'LineWidth',1.5,'XTick',1:4,'XTickLabel',xticks);
xlabel(xlab);
xlim([0.5 4.5]);
ylabel('Peak Firing Rate/Spontaneous Firing Rate - 1');
ylim([1.1*min(pdfr(:)) 1.1*max(pdfr(:))]);

legend([hb hf],{'ON Responses' 'OFF Responses' flab},'LineWidth',1.5,'Location','NorthWest');

%%

figName = sprintf('P%d_responses_first_%d_bins_no_double_dipping_separate_detection_n_spikes_spont_detection',P,k);
saveas(gcf,figName,'fig');
export_fig(figName,'-eps','-png','-transparent','-painters');

return;

%%

commonBestUnits = intersect(allBestUnits{3},allBestUnits{4});
n = numel(commonBestUnits);

efr = zeros(n,2);
sfr = zeros(n,2);

for jj = 1:2
    idx = ismember(allBestUnits{jj+2},commonBestUnits);
    efr(:,jj) = evokedFR{jj+2,1}(idx);
    sfr(:,jj) = spontFR{jj+2,1}(idx);
end