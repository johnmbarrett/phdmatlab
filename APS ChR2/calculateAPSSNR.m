allAPSExperiments = [47 49 50 54 56:62]; % exclude 63 because just generally a crap experiment
nAPSExperiments = numel(allAPSExperiments);

Ns = cell(nAPSExperiments,1);
SNRs = cell(nAPSExperiments,1);

currentDir = pwd;

for hh = 1:nAPSExperiments
    try
        exptDir = sprintf('JBOG%04d',allAPSExperiments(hh));
        cd(exptDir);
        
        ctrlData = dir('ChR2_cells_*_ctrl.mat');
        ctrlData = {ctrlData.name};

        drugData = dir('ChR2_cells_*_drug.mat');
        drugData = {drugData.name};

        nRecordings = numel(ctrlData);
        assert(numel(drugData) == nRecordings);

        allData = [ctrlData' drugData'];

        nCells = zeros(nRecordings,2);
        muSignal = cell(nRecordings,2);
        muNoise = cell(nRecordings,2);
        sigmaNoise = cell(nRecordings,2);

        for ii = 1:nRecordings
            for jj = 1:2
                load(allData{ii,jj},'allNSpikes','bestUnits','goodUnits','samples');

                nCells(ii,jj) = numel(bestUnits);
                muSignal{ii,jj} = zeros(nCells(ii,jj),1);
                muNoise{ii,jj} = zeros(nCells(ii,jj),1);
                sigmaNoise{ii,jj} = zeros(nCells(ii,jj),1);

                for kk = 1:nCells(ii,jj)
                    sample = samples{ismember(goodUnits,bestUnits(kk))};

                    muNoise{ii,jj}(kk) = mean(sample);
                    sigmaNoise{ii,jj}(kk) = std(sample);

                    muSignal{ii,jj}(kk) = max(mean(allNSpikes{kk}{1}));
                end
            end
        end

        Ns{hh} = nCells;
        SNRs{hh} = cellfun(@(ms,mn,sn) (ms-mn)./sn,muSignal,muNoise,sigmaNoise,'UniformOutput',false);
        
        cd(currentDir);
    catch err
        logMatlabError(err);
        cd(currentDir);
    end
end

%%

NDs = {                 ...
        [2.0 2.2 2.2];  ...
        [2.0 2.2 2.2];  ...
        [2.0 2.2 2.2];  ...
        [2.2 2.2 1.9];  ...
        [2.2 2.2 1.9];  ...
        [2.2 2.2 1.9];  ...
        [2.2 2.2 1.9];  ...
        [2.2 2.2 1.9];  ...
        [2.2 2.2];      ...
        [2.2 2.2 1.9];  ...
        [2.2 1.9 2.2];  ...
        [2.2 2.2 1.9]   ...
    };

n = sum(cellfun(@(ND) numel(unique(ND)),NDs));

NC = zeros(n,2);
ND = zeros(n,2);
SNR = zeros(n,2);

mSNR = cellfun(@(S) cellfun(@median,S),SNRs,'UniformOutput',false);

nn = 0;
for ii = 1:nAPSExperiments
    [uND,~,iND] = unique(NDs{ii});
    
    for jj = 1:max(iND)
        nn = nn + 1;
        
        ND(nn,:) = uND(jj);
        
        for kk = 1:2
            NC(nn,kk) = median(Ns{ii}(iND == jj,kk));
            SNR(nn,kk) = median(mSNR{ii}(iND == jj,kk));
        end
    end
end

%%

save('aps_snr_data.mat','Ns','SNRs','NDs','NC','ND','SNR');

%%

NC = reshape(NC,2*n,1);
ND = reshape(ND,2*n,1);
SNR = reshape(SNR,2*n,1);
D = kron([0;1],ones(n,1));

% xlswrite('aps_snr_data.xlsx',{'Num Cells' 'SNR' 'Neutral Density' 'Drug'},'A1:D1');
% xlswrite('aps_snr_data.xlsx',[NC SNR ND D],sprintf('A2:D%d',2*n+1));

%%

[uND,~,iND] = unique(ND);

means = zeros(max(iND),2,2);
stds = zeros(max(iND),2,2);

col = 0;
for ii = 1:max(iND)
    for jj = 0:1
        col = col + 1;
        
        colLetter = char('A'+col-1);
        
%         xlswrite('aps_snr_data',{sprintf('ND %d Drug %d',ii,jj)},'Num Cells',sprintf('%s1:%s1',colLetter,colLetter));
%         xlswrite('aps_snr_data',{sprintf('ND %d Drug %d',ii,jj)},'SNR',sprintf('%s1:%s1',colLetter,colLetter));
        
        data = NC(iND == ii & D == jj);
        means(ii,jj+1,1) = mean(data);
        stds(ii,jj+1,1) = std(data);
        
%         xlswrite('aps_snr_data',data,'Num Cells',sprintf('%s2:%s%d',colLetter,colLetter,numel(data)+1));
        
        data = SNR(iND == ii & D == jj);
        means(ii,jj+1,2) = mean(data);
        stds(ii,jj+1,2) = std(data);
        
%         xlswrite('aps_snr_data',data,'SNR',sprintf('%s2:%s%d',colLetter,colLetter,numel(data)+1));
    end
end

%%

figure;
ylabels = {'# Responsive Cells' 'Signal-to-Noise Ratio'};

for ii = 1:2
    subplot(1,2,ii);
    
    [hb,he] = barwitherr(flipud(stds(:,:,ii)),flipud(means(:,:,ii)));
    set([hb he],'LineWidth',1.5);
    set(gca,'LineWidth',1.5,'XTickLabel',[2.2 2.0 1.9]);
    xlabel('Neutral Density Filter');
    xlim([0.5 3.5])
    ylabel(ylabels{ii});
end

legend({'Control' 'MFA'},'LineWidth',1.5,'Location','NorthWest');

%%

saveas(gcf,'aps_snr_data','fig');
saveas(gcf,'aps_snr_data','emf');
export_fig('aps_snr_data','-eps','-png',-'m4','-transparent','-painters');