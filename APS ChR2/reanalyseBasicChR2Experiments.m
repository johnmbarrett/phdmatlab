% P = 95;
P = 96;

%%

load('ChR2_cells');
load('EventNo');
load('spiketimestamps.mat');

if P == 96
    EventNo = EventNo([1 2 3 5]);
end

%%

n = numel(bestUnits);
spontFR = zeros(n,4);

for ii = 1:4
    if ii == 1
        t0 = 0;
    else
        t0 = EventNo{ii-1}(end);
    end
    
    t1 = EventNo{ii}(1);
    T = t1-t0;
    
    spontFR(:,ii) = cellfun(@(t) sum(t > t0 & t < t1),spiketimestamps(bestUnits))/T;
end

bad = min(spontFR,[],2) == 0;
disp(n);
disp(sum(bad));
disp(n-sum(bad));

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

peakDFR = peakFR./repmat(spontFR,[1 1 2]);
peakDFR = peakDFR(~bad,:,:);

%%

pdfr = prctile(peakDFR,[25 50 75]);
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
    yfill = [0.1 100 100 0.1];
    ylims = [0 70];
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
ylabel('Peak Firing Rate/Spontaneous Firing Rate');
ylim(ylims);

legend([hb hf],{'ON Responses' 'OFF Responses' flab},'LineWidth',1.5,'Location','NorthWest');

%%

figName = sprintf('P%d_responses',P);
saveas(gcf,figName,'fig');
export_fig(figName,'-eps','-png','-transparent','-painters');