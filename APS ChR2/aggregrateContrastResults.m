contrastExperiments = [64 65 67:69];
nContrastExperiments = numel(contrastExperiments);

preStimss = zeros(size(contrastExperiments));
preStimss(end) = 2;

topDir = pwd;

skipDetection = true;
reanalyse = false;
lumstepAnalysisSuffix = '_inc_inf';
constepAnalysisSuffix = '';
revgratAnalysisSuffix = '_250ms_before_after_abs_diff';

%%

if reanalyse
    for ii = 1:nContrastExperiments
        exptDir = sprintf('JBOG%04d',contrastExperiments(ii));
        cd(exptDir);

        try
            preStims = preStimss(ii);

            analyseAPSChR2ContrastExperiment;

            cd(topDir);

        catch err
            logMatlabError(err);
            cd(topDir);
        end
    end
end

%%

allLumDSpikes = zeros(nContrastExperiments,8,2);
allLumPSpikes = zeros(nContrastExperiments,8,2);
allLumSNR = zeros(nContrastExperiments,8,2);
allConDSpikes = zeros(nContrastExperiments,25,2);
allF2F1Ret = zeros(nContrastExperiments,6,2);

for ii = 1:nContrastExperiments
    exptDir = sprintf('JBOG%04d',contrastExperiments(ii));
    
    lumstepFile = sprintf('%s\\luminance_steps_responses%s.mat',exptDir,lumstepAnalysisSuffix);
    
    L = load(lumstepFile,'nSpikes','pSpikes','snr','spontFR');
    
    constepFile = sprintf('%s\\contrast_steps_responses%s.mat',exptDir,constepAnalysisSuffix);
    
    C = load(constepFile,'dSpikes');
    
    for jj = 1:2
        allLumDSpikes(ii,:,jj) = median(4*L.nSpikes{jj}-repmat(L.spontFR{jj}',8,1),2);
        allLumPSpikes(ii,:,jj) = median(L.pSpikes{jj},2);
        allLumSNR(ii,:,jj) = median(L.snr{jj},2);
%         allConDSpikes(ii,:,jj) = median(C.dSpikes{jj},2);
%         allLumPSpikes(ii,:,jj) = median(pSpikes{jj}(:,all(isfinite(pSpikes{jj}),1)),2);
%         allLumSNR(ii,:,jj) = median(snr{jj}(:,all(isfinite(snr{jj}),1),:),2);
        allConDSpikes(ii,:,jj) = median(C.dSpikes{jj}(:,all(isfinite(C.dSpikes{jj}),1)),2);
    end
    
    revgratFile = sprintf('%s\\reversing_gratings_responses%s.mat',exptDir,revgratAnalysisSuffix);
    
    R = load(revgratFile,'f2f1Retina');
    
    allF2F1Ret(ii,:,:) = R.f2f1Retina;
end

%%

conditions = 0; % suppress weird error
load(sprintf('%s\\constep sequence.mat',exptDir));

initLums = lums(conditions(2:end-1,2));
contrasts = 100*(lums(conditions(2:end-1,3))./initLums-1);
[~,sortIndices] = sort(contrasts);
% initLums = initLums(sortIndices);
[~,~,il] = unique(initLums);
[~,~,ic] = unique(contrasts);

load(sprintf('%s\\lumstep sequence.mat',exptDir));

%%

figure;
% medianErrorbar(lums,allLumPSpikes);
errorbar(repmat(lums(:),1,2),squeeze(mean(allLumDSpikes)),squeeze(std(allLumDSpikes)));

%%

figure;
% medianErrorbar(lums,allLumPSpikes);
errorbar(repmat(lums(:),1,2),squeeze(mean(allLumPSpikes)),squeeze(std(allLumPSpikes)));

%%

figure;
% medianErrorbar(lums,allLumSNR);
errorbar(repmat(lums(:),1,2),squeeze(mean(allLumSNR)),squeeze(std(allLumSNR)));

%%

figure;
% medianErrorbar(contrasts,allConDSpikes(:,sortIndices,:));
errorbar(repmat(contrasts(sortIndices)',1,2),squeeze(mean(allConDSpikes(:,sortIndices,:))),squeeze(std(allConDSpikes(:,sortIndices,:))));
figFile = ['aps_constep_results' constepAnalysisSuffix];
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');

%%

load(sprintf('%s\\revgrat sequence.mat',exptDir));

%%

figure;
medianErrorbar(repmat((0.1:0.1:0.6)',1,2),allF2F1Ret(all(all(isfinite(allF2F1Ret),2),3),:,:));

legend({'Control' 'MFA'});
set(gca,'XTick',0.1:0.1:0.6,'XTickLabel',arrayfun(@(c) sprintf('%02.0f',100*c),C,'UniformOutput',false));
xlabel('Michelson Contrast (%)');
xlim([0 0.7]);
ylabel('|log(Reversal Response/Grating Baseline)|');
% ylim([0 1]);

figFile = ['aps_revgrat_results' revgratAnalysisSuffix];
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');

%%

save(['aps_contrast_data' constepAnalysisSuffix],'allLumDSpikes','allLumPSpikes','allLumSNR','allConDSpikes','initLums','contrasts','lums');

%%

contrastResultsFile = ['apschr2_spss_contrast_data' constepAnalysisSuffix '.xlsx'];
xlswrite(contrastResultsFile,reshape(permute(allLumPSpikes,[1 3 2]),nContrastExperiments,16),'% Change in Firing (Lum Steps)',sprintf('A1:P%d',nContrastExperiments));
xlswrite(contrastResultsFile,reshape(permute(allLumSNR,[1 3 2]),nContrastExperiments,16),'SNR (Luminance Steps)',sprintf('A1:P%d',nContrastExperiments));

%%

% xlswrite(contrastResultsFile,nan(6,212),'% Change in FR (Con Steps) 2','A1:HE6');
% xlswrite(contrastResultsFile,ones(6,1),'% Change in FR (Con Steps) 2','A1:A6');
% 
% for ii = 1:25
%     for jj = 1:2
%         col = 42*(il(ii)-1)+2*(ic(ii)-1)+jj+1;
%         excelCol = getExcelColumn(col);
%         xlswrite(contrastResultsFile,{sprintf('L%dC%dD%d',il(ii),ic(ii),jj)},'% Change in FR (Con Steps) 2',sprintf('%s1:%s1',excelCol,excelCol));
%         xlswrite(contrastResultsFile,allConDSpikes(:,ii,jj),'% Change in FR (Con Steps) 2',sprintf('%s2:%s%d',excelCol,excelCol,nContrastExperiments+1));
%     end
% end
% 
% xlswrite(contrastResultsFile,ones(6,1),'% Change in FR (Con Steps) 2','HE1:HE6');

%%

xlswrite(contrastResultsFile,reshape(permute(allConDSpikes,[1 3 2]),nContrastExperiments,50),'dFR',sprintf('A1:AX%d',nContrastExperiments));
xlswrite(contrastResultsFile,reshape(contrasts,5,5),'dFR',sprintf('AY1:BC%d',nContrastExperiments));

%%

spssData = zeros(250,5);
spssData(:,1) = repmat((1:5)',50,1);
spssData(:,2) = repmat(kron(initLums',ones(5,1)),2,1);
spssData(:,3) = repmat(kron(contrasts',ones(5,1)),2,1);
spssData(:,4) = kron((1:2)',ones(125,1));
spssData(:,5) = allConDSpikes(:);
spssData(:,6) = repmat(kron((1:5)',ones(5,1)),10,1);