% movingBarExperiments = 47:50;
% nMovingBarExperiments = numel(movingBarExperiments);
% 
% barPerf = zeros(2,2,nMovingBarExperiments);
% 
% for ii = 1:nMovingBarExperiments
%     exptDir = sprintf('JBOG%04d',movingBarExperiments(ii));
%     
%     load([exptDir '/bar_responses.mat']);
%     
%     barPerf(:,:,ii) = decoderPerformance;
% end
% 
% load([exptDir '/bars sequence.mat']);

% 48 is the one with the low potassium in the drug condition
oldGratingExperiments = [47 49 50 54 56];
nOldGratingExperiments = numel(oldGratingExperiments);

oldFreqPerf = zeros(6,2,nOldGratingExperiments);
oldContPerf = zeros(6,2,nOldGratingExperiments);

for ii = 1:nOldGratingExperiments
    exptDir = sprintf('JBOG%04d',oldGratingExperiments(ii));
    
    load([exptDir '/frequency_responses.mat']);
    
    oldFreqPerf(:,:,ii) = decoderPerformance;
    
    load([exptDir '/contrast_responses.mat']);
    
    oldContPerf(:,:,ii) = decoderPerformance;
end

load([exptDir '/frequency sequence.mat']);
load([exptDir '/contrast sequence.mat']);

%%

freqs = (30./(barWidths*4*2));

%%

figure;

errorbar(repmat(freqs',1,2),mean(100*oldFreqPerf,3),std(100*oldFreqPerf,[],3),'LineWidth',1.5);
line([0.01 1],[12.5 12.5],'Color','k','LineStyle','--','LineWidth',1.5);
legend({'Control' 'Drug' 'Chance'});
set(gca,'LineWidth',1.5,'XScale','log','XTick',wrev(freqs),'XTickLabel',arrayfun(@(x) sprintf('%0.3f',x),wrev(freqs),'UniformOutput',false));
xlabel('Spatial Frequency (cpd)');
xlim([10^-2 10^-0.4]);
ylabel('Decoder Performance (%)');
ylim([0 100]);

%%

saveas(gcf,'all_retinas_old_freq_perf','fig');
export_fig('all_retinas_old_freq_perf','-eps','-png','-transparent','-painters','-m4');
close(gcf);

%%

load('V:\retina\John B\phd backup\matlab\stimuli\aps example\gamma.mat','q');
L = polyval(q,0:255);
contrasts = (L(V2+1)-L(V1+1))./(L(V2+1)+L(V1+1));

%%

figure;
hold on;
errorbar(repmat((4:6)',1,2),mean(100*oldContPerf(1:3,:,:),3),std(100*oldContPerf(1:3,:,:),[],3),'LineWidth',1.5);
errorbar(repmat([3;2;1],1,2),mean(100*oldContPerf(4:6,:,:),3),std(100*oldContPerf(4:6,:,:),[],3),'LineWidth',1.5);
line([0 7],[12.5 12.5],'Color','k','LineStyle','--','LineWidth',1.5);
legend({'Control' 'Drug' 'Chance'});
set(gca,'LineWidth',1.5,'XTick',1:6,'XTickLabel',arrayfun(@(x) sprintf('%2.1f',x),100*contrasts([6 5 4 1 2 3]),'UniformOutput',false));
xlabel('Michelson Contrast (%)');
xlim([0.5 6.5]);
ylabel('Decoder Performance (%)');
ylim([0 100]);

%%

X = zeros(6*nOldGratingExperiments,2);
datas = {oldFreqPerf oldContPerf};

ps = zeros(1,2);
tables = cell(1,2);
statss = cell(1,2);

for ii = 1:2
    for jj = 1:6
        for kk = 1:nOldGratingExperiments
            X(nOldGratingExperiments*(jj-1)+kk,:) = datas{ii}(jj,:,kk);
        end
    end
    
    [ps(ii),tables{ii},statss{ii}] = friedman(X,nOldGratingExperiments,'off');
end

%%

Y = zeros(nOldGratingExperiments,12);
infixes = {'freq' 'cont'};

for hh = 1:2
    for ii = 1:6
        for jj = 1:2
            Y(:,2*(ii-1)+jj) = datas{hh}(ii,jj,:);
        end
    end
    
    xlswrite('apschr2_spss_data.xlsx',Y,sprintf('old %s perf',infixes{hh}),sprintf('A1:L%d',nOldGratingExperiments));
end

%%

saveas(gcf,'all_retinas_old_cont_perf','fig');
export_fig('all_retinas_old_cont_perf','-eps','-png','-transparent','-painters','-m4');
close(gcf);

%%

freqGratingExperiments = 57:63; % 63 was crap
contGratingExperiments = 57:61; % 62 is where I started on the luminance gratings
nFreqGratingExperiments = numel(freqGratingExperiments);
nContGratingExperiments = numel(contGratingExperiments);

freqPerf = zeros(4,6,2,nFreqGratingExperiments);
contPerf = zeros(4,6,2,nContGratingExperiments);

data = zeros(24,3,2);

for ii = 1:nFreqGratingExperiments
    exptDir = sprintf('JBOG%04d',freqGratingExperiments(ii));
    
    load([exptDir '/frequency_2afc_responses.mat']);
    load([exptDir '/frequency sequence.mat']);
    
    freqs = (30./(barWidths*4*2));
    data(:,1,:) = repmat(-log(freqs'),[4 1 2]);
    data(:,3,:) = 50;
    
    for jj = 1:6
        for kk = 1:4
            for ll = 1:2
                data(6*(kk-1)+jj,2,ll) = sum(detectionSuccess(:,kk,jj,ll));
            end
        end
    end
    
    psignifitFile = sprintf('JBOG%04d_freq_psignifit_data.txt',freqGratingExperiments(ii));
    fout = fopen(psignifitFile,'w');
    
    for jj = 1:24
        fprintf(fout,'%d\t%d\t%d\t%d\r\n',data(jj,1,1),data(jj,2,1),data(jj,2,2),data(jj,3,1));
    end
    
    fclose(fout);
    
    freqPerf(:,:,:,ii) = detectionPerf;
end

%%
    
allData = zeros(24,3,2);
allData(:,1,1) = repmat(-log(freqs'),4,1);
allData(:,1,2) = repmat(-log(freqs'),4,1);
allData(:,3,:) = 50;

for hh = 1 %:6
    for jj = 1:6
        for kk = 1:4
            for ll = 1:2
                allData(24*(hh-1)+6*(kk-1)+jj,2,ll) = mean(freqPerf(kk,jj,ll,hh),4)/100;
            end
        end
    end
end

psignifitFile = 'all_retinas_freq_psignifit_data.txt';
fout = fopen(psignifitFile,'w');

for jj = 1:24
    fprintf(fout,'%d\t%d\t%d\t%d\r\n',allData(jj,1,1),allData(jj,2,1),allData(jj,2,2),allData(jj,3,1));
end

fclose(fout);

%%

for ii = 1:nContGratingExperiments
    exptDir = sprintf('JBOG%04d',contGratingExperiments(ii));
    
    load([exptDir '/contrast_2afc_responses.mat']);
    
    contPerf(:,:,:,ii) = detectionPerf;
end

%%


load([exptDir '/contrast sequence.mat']);

%%

P = repmat(((0:3)/4)',[1 6 2 nFreqGratingExperiments]);

F = log(repmat(freqs,[4 1 2 nFreqGratingExperiments]));

D = repmat(reshape([0 1],[1 1 2 1]),[4 6 1 nFreqGratingExperiments]);

R = repmat(reshape(1:nFreqGratingExperiments,[1 1 1 nFreqGratingExperiments]),[4 6 2 1]);

Y = freqPerf;

%%

b = glmfit([P(:) F(:) R(:)],Y(:),'normal');
yfit = glmval(b,[P(:) F(:) R(:)],'identity');
figure;
[fsort,si] = sort(F(:));
plot(F(:),Y(:),'o',fsort,yfit(si),'-');

%%

phaseAverageFreqPerf = squeeze(mean(freqPerf,1));

%%

X = zeros(6*nFreqGratingExperiments,2);

for ii = 1:6
    for jj = 1:nFreqGratingExperiments
        X(nFreqGratingExperiments*(ii-1)+jj,:) = phaseAverageFreqPerf(ii,:,jj);
    end
end

[pf,tablef,statsf] = friedman(X,nFreqGratingExperiments);

%%

Y = zeros(nFreqGratingExperiments,12);

for ii = 1:6
    for jj = 1:2
        Y(:,2*(ii-1)+jj) = phaseAverageFreqPerf(ii,jj,:)/100;
    end
end
    
xlswrite('apschr2_spss_data.xlsx',Y,'new freq perf',sprintf('A1:L%d',nFreqGratingExperiments));

%%

figure;

errorbar(repmat(freqs',1,2),mean(phaseAverageFreqPerf,3),std(phaseAverageFreqPerf,[],3),'LineWidth',1.5);
line([0.01 1],[50 50],'Color','k','LineStyle','--','LineWidth',1.5);
legend({'Control' 'Drug' 'Chance'});
set(gca,'LineWidth',1.5,'XScale','log','XTick',wrev(freqs),'XTickLabel',arrayfun(@(x) sprintf('%0.3f',x),wrev(freqs),'UniformOutput',false));
xlabel('Spatial Frequency (cpd)');
xlim([10^-2 10^-0.4]);
ylabel('Decoder Performance (%)');
ylim([40 100]);

%%

saveas(gcf,'all_retinas_freq_perf','fig');
export_fig('all_retinas_freq_perf','-eps','-png','-transparent','-painters','-m4');
close(gcf);

%%

X = log(repmat(freqs,4,1));

freqParams = zeros(2,2,nFreqGratingExperiments);

colors = 'bg';

sigmoid = @(p,x) 1-0.5./(1+exp(-p(1)*(x-p(2))));

% figure;
% 
% fplot(@(x) sigmoid([1 -3],x),[-5 0])

%%

for hh = 1:nFreqGratingExperiments
    figure;
    
    for ii = 1:2
        Y = squeeze(freqPerf(:,:,ii,hh))/100;
        freqParams(:,ii,hh) = lsqcurvefit(sigmoid,[1 -3],X(:),Y(:),[0 -Inf],[Inf Inf]);
        subplot(1,2,ii);
        hold on;
        plot(X(:),Y(:),'Color',colors(ii),'LineStyle','none','Marker','o');
        fplot(@(x) sigmoid(freqParams(:,ii,hh),x),xlim,'Color',colors(ii));
    end
end

%%

ydatas = permute(mean(mean(contPerf,4),1),[2 3 1])/100;
% ydatas = permute(mean(contPerf,4),[2 1 3])/100;
% ydatas = reshape(ydatas(:,:,2),24,1);

nonlcon = @contrastGratingDecoderSimulationNonlinearConstraints;
objFun = @(p) sum((ydatas(:,2)-mean(mean(reshape(cell2mat(arrayfun(@(n) contrastGratingDecoderSimulation(p(1),p(2),p(3),p(4)),1:10,'UniformOutput',false)),[6 4 10]),3),2)).^2);
options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','Aineq',[1 1 0 0],'bineq',80,'lb',[0 0 0 0],'ub',[Inf 1 Inf Inf],'nonlcon',nonlcon,'objective',objFun,'x0',[8/3 -5/3 1 2.5],'options',options);

[params,fval,exitflat,output,solutions] = fmincon(objFun,[1 0 1 2.5],[1 1 0 0],80,[],[],[0 0 0 0],[Inf 1 Inf Inf],nonlcon,options);

figure;
plot(0.1:0.1:0.6,ydatas(:,2),0.1:0.1:0.6,mean(contrastGratingDecoderSimulation(params(1),params(2),params(3),params(4)),2));

% gs = GlobalSearch;
% [params,fval,exitflat,output,solutions] = run(gs,problem);

%%

phaseAverageContPerf = squeeze(mean(contPerf,1));

%%

X = zeros(6*nContGratingExperiments,2);

for ii = 1:6
    for jj = 1:nContGratingExperiments
        X(nContGratingExperiments*(ii-1)+jj,:) = phaseAverageContPerf(ii,:,jj);
    end
end

[pc,tablec,statsc] = friedman(X,nContGratingExperiments);

%%

Y = zeros(nContGratingExperiments,12);

for ii = 1:6
    for jj = 1:2
        Y(:,2*(ii-1)+jj) = phaseAverageContPerf(ii,jj,:)/100;
    end
end
    
xlswrite('apschr2_spss_data.xlsx',Y,'new cont perf',sprintf('A1:L%d',nContGratingExperiments));

%%

figure;
conts = 10:10:60;
errorbar(repmat(conts',1,2),mean(phaseAverageContPerf,3),std(phaseAverageContPerf,[],3),'LineWidth',1.5);
line([0 70],[50 50],'Color','k','LineStyle','--','LineWidth',1.5);
legend({'Control' 'Drug' 'Chance'});
set(gca,'LineWidth',1.5,'XTick',conts,'XTickLabel',conts);
xlabel('Michelson Contrast (%)');
xlim([0 70]);
ylabel('Decoder Performance (%)');
ylim([30 100]);

%%

saveas(gcf,'all_retinas_cont_perf','fig');
export_fig('all_retinas_cont_perf','-eps','-png','-transparent','-painters','-m4');
close(gcf);

%%

X = repmat(conts,4,1);

contParams = zeros(2,2,nFreqGratingExperiments);

colors = 'bg';

sigmoid = @(p,x) 0.5./(1+exp(-p(1)*(x-p(2))))+0.5;

% figure;
% 
% fplot(@(x) sigmoid([1 -3],x),[-5 0])

%%

for hh = 1:nFreqGratingExperiments
    figure;
    
    for ii = 1:2
        Y = squeeze(contPerf(:,:,ii,hh))/100;
        contParams(:,ii,hh) = lsqcurvefit(sigmoid,[1 35],X(:),Y(:),[0 -Inf],[Inf Inf]);
        subplot(1,2,ii);
        hold on;
        plot(X(:),Y(:),'Color',colors(ii),'LineStyle','none','Marker','o');
        fplot(@(x) sigmoid(contParams(:,ii,hh),x),xlim,'Color',colors(ii));
    end
end

%%

letterExperiments = [56:59 61 62 63];
nLetterExperiments = numel(letterExperiments);

posPerf = zeros(15,2,nLetterExperiments);
scalePerf = cell(4,2,nLetterExperiments);

data = zeros(15,3,2);
data(:,3,:) = 25;

for ii = 1:nLetterExperiments
    exptDir = sprintf('JBOG%04d',letterExperiments(ii));
    
    load([exptDir '/letter_responses.mat']);
    
    posPerf(:,:,ii) = perfPerPosition;
    
    load([exptDir '/letters sequence.mat']);
    
    cumPositionsPerScale = [0 cumsum(positionsPerScale)];
    
    for jj = 1:nScales
        idx = cumPositionsPerScale(jj)+1:cumPositionsPerScale(jj+1);
        for kk = 1:2
            scalePerf{jj,kk,ii} = perfPerPosition(idx,kk);
            data(idx,1,kk) = 4*(sizes(jj)/(30*5));
            data(idx,2,kk) = perfPerPosition(idx,kk);
        end
    end
    
    outputFile = sprintf('JBOG%04d_letters_psignifit_data.txt',letterExperiments(ii));
    fout = fopen(outputFile,'w');
    
    for jj = 1:15
        fprintf(fout,'%d\t%d\t%d\t%d\r\n',data(jj,1,1),data(jj,2,1),data(jj,2,2),data(jj,3,1));
    end
       
    fclose(fout);
end

%%

% psignifitResults = importdata('letters_psignifit_results_small_lambda_gauss_m_prior.txt');
psignifitResults = importdata('letters_psignifit_results_free_lambda_gauss_16_8_m_prior.txt');
psignifitResults = psignifitResults.data;

midpoints = psignifitResults(:,[2 6]);
widths = psignifitResults(:,[3 7]);
lambdas = psignifitResults(:,[4 8]);

%%

figure;
hold on;

colours = 'bg';
sigmoidFun = @(m,w,l,a) @(x) 0.1+(0.9-l)./(1+exp(-2.*log(1/a-1).*(x-m)./w));
inverseFun = @(m,w,l,a) @(y) m-w.*log((0.9-l)./(y-0.1)-1)./(2*log(1./a-1));

maxLogMar = ceil(10*log(160*60/5)/log(10))/10;

thresh20s = zeros(nLetterExperiments,2);
logMARs = maxLogMar*ones(nLetterExperiments,2)+0.1;

for ii = 1:2
    for jj = 1:nLetterExperiments
        sigmoid = sigmoidFun(midpoints(jj,ii),widths(jj,ii),lambdas(jj,ii),0.1);
        fplot(sigmoid,[0 20],colours(ii));
        inverse = inverseFun(midpoints(jj,ii),widths(jj,ii),lambdas(jj,ii),0.1);
        thresh20s(jj,ii) = ceil(10*log(inverse(0.2)*60)/log(10))/10;
        
        logMarRange = thresh20s(jj,ii):0.1:maxLogMar;
        
        if isempty(logMarRange)
            continue;
        end
        
        angleRange = (10.^logMarRange)/60;
        logMARs(jj,ii) = maxLogMar+0.1-0.02*sum(5*sigmoid(angleRange));
    end
end

ylim([0 1])

snellenDenominators = 20*10.^logMARs;

%%

positionAverageScalePerf = 100*squeeze(cellfun(@mean,scalePerf));

%%

X = zeros(6*nLetterExperiments,2);

for ii = 1:4
    for jj = 1:nLetterExperiments
        X(nLetterExperiments*(ii-1)+jj,:) = positionAverageScalePerf(ii,:,jj);
    end
end

[pl,tablel,statsl] = friedman(X,nLetterExperiments);

%%

Y = zeros(nLetterExperiments,8);

for ii = 1:4
    for jj = 1:2
        Y(:,2*(ii-1)+jj) = positionAverageScalePerf(ii,jj,:)/100;
    end
end
    
xlswrite('apschr2_spss_data.xlsx',Y,'letters perf',sprintf('A1:L%d',nLetterExperiments));

%%

letterSizes = sizes*4/(30*5);

%%

figure;
errorbar(repmat(letterSizes',1,2),mean(positionAverageScalePerf,3),std(positionAverageScalePerf,[],3),'LineWidth',1.5);
line([0 20],[10 10],'Color','k','LineStyle','--','LineWidth',1.5);
legend({'Control' 'Drug' 'Chance'},'Location','NorthWest');
set(gca,'LineWidth',1.5,'XTick',wrev(letterSizes),'XTickLabel',arrayfun(@(x) sprintf('%4.2f',x),wrev(letterSizes),'UniformOutput',false));
xlabel('Feature Size (degrees of visual angle)');
xlim([2 17]);
ylabel('Decoder Performance (%)');
ylim([0 100]);

%%

save('aps_acuity_data.mat','freqs','freqPerf','letterSizes','posPerf','scalePerf');

%%

saveas(gcf,'all_retinas_letter_perf','fig');
export_fig('all_retinas_letter_perf','-eps','-png','-transparent','-painters','-m4');
close(gcf);

%%

X = cell2mat(arrayfun(@(l,n) repmat(l,1,n),letterSizes,positionsPerScale,'UniformOutput',false))';

letterParams = zeros(2,2,nFreqGratingExperiments);

colors = 'bg';

sigmoid = @(p,x) 0.9./(1+exp(-p(1)*(x-p(2))))+0.1;

% figure;
% 
% fplot(@(x) sigmoid([1 -3],x),[-5 0])

%%

for hh = 1:nLetterExperiments
    figure;
    
    for ii = 1:2
        Y = posPerf(:,ii,hh);
        letterParams(:,ii,hh) = lsqcurvefit(sigmoid,[1 10],X(:),Y(:),[0 0],[Inf Inf]);
        subplot(1,2,ii);
        hold on;
        plot(X(:),Y(:),'Color',colors(ii),'LineStyle','none','Marker','o');
        fplot(@(x) sigmoid(letterParams(:,ii,hh),x),xlim,'Color',colors(ii));
    end
end

%%

save('H:\apschr2_data.mat','freqs','oldFreqPerf','contrasts','oldContPerf','freqPerf','conts','contPerf','sizes','scalePerf');