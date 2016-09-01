%% INITIALISATION

[numdata,txtdata] = xlsread('responsive cell recordings.xlsx','A2:E22');

recordings = numdata(:,1);
middleVolts = numdata(:,4);
offsets = numdata(:,5);

drugs = txtdata(:,1);

uniqueDrugs = unique(drugs);

retinas = txtdata(:,2);

%% REDO ANALYSIS

nRecordings = numel(recordings);
disp(nRecordings);
topDir = pwd;
concs = [0 0 0 0 0 0 10 20 40 80 80 80 80 80 80 0]';

% moving bar analysis only : 3 5 8 13 14
for ii = 1:nRecordings
    middleVolt = middleVolts(ii);
    
    volts = middleVolt*ones(16,1);
    
    for jj = 0:1
        volts((1:5)+jj*9) = (6:10)';
    end
    
    detectionRecordings = middleVolt+[-5 3];
    
    recording = recordings(ii);
    recDir = sprintf('JBOG%04d',recording);
    cd(recDir);
    
    try
%         reanalyseExperiment(retinas{ii},drugs{ii},concs,volts,detectionRecordings,offsets(ii));
%         movingBarAnalysis(retinas{ii},detectionRecordings);
        uledSTA(retinas{ii},'flashing squares',detectionRecordings);
%         movingBarAnalysisTemp(retinas{ii},detectionRecordings,1,true,false,'mbg_allresp_fullprior',true)
        close all;
        cd(topDir);
    catch err
        logMatlabError(err);
        cd(topDir);
    end
end

assert(false);

%% AGGREGATE STIMULATION DATA

aggregateStimData(recordings,retinas,drugs,middleVolts);

%% LOAD AGGREGATED STIM DATA

load('all_stim_data_mbg_all_resp_full_prior.mat');

%% LOAD AGGREGATED SPONT DATA

load('og_aggregate_data.mat');

spontPs = zeros(3,2);
data = permute(cat(4,data{:}),[1 2 4 3]);

for ii = 1:2
    for jj = 1:3
        X = data(:,:,jj,ii);
        X = X(isfinite(X(:,1)),:);
        spontPs(jj,ii) = friedman(X,1,'off');
    end
end

%% STATS
allPs = [spontPs(:,1); spontPs(:,2); nPs; snrPs; threshPs; barPs];
% allPs = [spontPs(:,1); spontPs(:,2); nPs; snrPs; threshPs; sdPs; eccentricityPs; barPs];
[barhs,~,as] = holmBonferroni(allPs);

labels = cellfun(@(l) cellfun(@(d) sprintf('% 5s vs %s',d,l),{'18BGA' 'Flu' 'MFA'}','UniformOutput',false),{'Spont Firing' 'Spont Oscillations' 'Responsive Cells' 'SNR' 'Thresholds' 'Bar Performance'}','UniformOutput',false);
% labels = cellfun(@(l) cellfun(@(d) sprintf('% 5s vs %s',d,l),{'18BGA' 'Flu' 'MFA'}','UniformOutput',false),{'Spont Firing' 'Spont Oscillations' 'Responsive Cells' 'SNR' 'Thresholds' 'RF SD' 'RF Eccentricity' 'Bar Performance'}','UniformOutput',false);
labels = vertcat(labels{:});

for ii = 1:numel(allPs)
    fprintf('%27s\t\t%f\t%f\t%d\n',labels{ii},allPs(ii),as(ii),barhs(ii));
end

%% PLOT SPONT DATA

ylims = [0 120; 0 180];
ylabels = {'Spontaneous Firing Rate' 'Oscillation Strength'};

figure;
set(gcf,'Position',[100 100 780 400]); %'PaperPositionMode','manual','PaperUnits','centimeters','PaperSize',[18 10],'PaperPosition',[0 0 18 10]);

for ii = 1:2
    subplot('Position',[0.075+0.5*(ii-1) 0.1 0.4 0.85]);
    pspont = prctile(100*data(:,:,:,ii),[25 50 75]);
    mspont = squeeze(pspont(2,:,:));
    lspont = mspont-squeeze(pspont(1,:,:));
    uspont = squeeze(pspont(3,:,:))-mspont;
    [barhs,errhs] = barwitherr(cat(3,lspont,uspont),mspont,'LineWidth',1.5);
    
    for jj = 1:numel(errhs)
        set(errhs(jj),'LineWidth',1.5);
    end
    
    if ii == 2
        legend(barhs,{'18BGA' 'Flu' 'MFA'},'FontSize',8,'Location','NorthWest');
    end
    
    set(gca,'FontSize',10,'LineWidth',1.5,'XTickLabel',{'Ctrl' '10' '20' '40' '80' 'Wash'});
    
    text(-0.5,ylims(ii,2),char('A'+ii-1),'FontSize',24,'FontWeight','bold');
    
    xlabel('Drug Concentration ({\mu}M)','FontSize',10);
    xlim([0.5 6.5])
    ylabel(sprintf('%s (%% of Control)',ylabels{ii}),'FontSize',10);
    ylim(ylims(ii,:));
end

%% SAVE SPONT DATA FIG FOR PRES

figPath = 'V:\\retina\\John B\\phd backup\\presentations\\seattle jun 15\\fig2_spont_data';
    
export_fig(figPath,'-m4','-png','-transparent','-painters');

%% SAVE SPONT DATA FIG FOR PUB

figPath = 'V:\\retina\\John B\\phd backup\\papers\\Frontiers Degeneration of Visual Processing CFP\\img\\fig2_spont_data';
    
export_fig(figPath,'-eps','-png','-transparent','-painters');
% print(gcf,figPath,'-depsc','-r600');
% print(gcf,figPath,'-dpng','-r600');
saveas(gcf,figPath,'fig');
close(gcf);

%% INDEXES FOR SUBANALYSES

nidx = cellfun(@(drug) 1:sum(strcmp(drug,drugs)),uniqueDrugs,'UniformOutput',false);
aidx = cellfun(@(drug) find(strcmp('A',retinas(strcmp(drug,drugs)))),uniqueDrugs,'UniformOutput',false);
bidx = cellfun(@(drug) find(strcmp('B',retinas(strcmp(drug,drugs)))),uniqueDrugs,'UniformOutput',false);
nsidxs = [nidx aidx bidx];
nsidxs{1,4} = nsidxs{1,1}(2:end-1);
snridxs = nsidxs;

for ii = 1:4
    snridxs{1,ii}(end) = []; % exclude last 18BGA
    snridxs{2,ii}(snridxs{2,ii} == 3) = []; % exclude bad Flu
    snridxs{2,ii}(snridxs{2,ii} > 3) = snridxs{2,ii}(snridxs{2,ii} > 3) - 1;
end

%% DATA FOR SUBANALYSIS PLOTS

ydatas = {ns threshData cellfun(@(X) permute(X(:,end,:),[1 3 2]),snrFinData,'UniformOutput',false)};
colours = 'rgb';
figNums = [3 4 6];
figSuffixes = {'ncells' 'thresh' 'snr'};
ylabels = {'# Responsive Cells (% of Control)' 'Threshold (ms)' 'Signal-to-Noise Ratio'};
titles = {'All Retinas' 'A Retinas' 'B Retinas' 'Only 8V'};

%% PLOT SUBANALYSIS FIGURES

subfig = figure;
pubfig = figure;
set(gcf,'Position',[100 100 368 246]);

for hh = 1 %:4
    for ii = 1:3
        if hh == 1
            figure(pubfig);
            clf;
            hold on;
        else
            subplot(3,3,(ii-1)*3+hh-1);
            hold on;
        end

        if ii == 1
            ydata = cellfun(@(X) 100*X./repmat(X(:,1),1,size(X,2)),ydatas{ii},'UniformOutput',false);
        else
            ydata = ydatas{ii};
        end

        barhs = zeros(3,1);

        for jj = 1:3
            if ii == 1
                idx = nsidxs{jj,hh};
            else
                idx = snridxs{jj,hh};
            end
            
            [barhs(jj),m] = medianErrorbar(1:6,ydata{jj}(idx,:),'Color',colours(jj),'LineWidth',1.5);

%             if ii == 3 && jj > 1
%                 line([4 5; 5 6],[m(4) 1e4; 1e4 m(6)],'Color',colours(jj),'LineStyle','-','LineWidth',1.5);
%             end
        end

        legend(barhs,{'18BGA' 'Flu' 'MFA'},'FontSize',10,'Location','NorthWest');
        set(gca,'FontSize',10,'LineWidth',1.5,'Position',[0.175 0.175 0.775 0.775],'XTick',1:6,'XTickLabel',{'Ctrl' '10' '20' '40' '80' 'Wash'});
        xlabel('Drug Concentration ({\mu}M)','FontSize',10);
        xlim([0.5 6.5]);
%         ylab = sprintf('%s (%% of control)',ylabels{ii});
        
%         if ii < 3
%             ylabel(ylab,'FontSize',10);
%         else
            ylabel(ylabels{ii});
%             % ylabel cut off slightly in PDF if use default ylabel settings
%             % for SNR
%             text(-0.5,550,ylab,'FontSize',10,'HorizontalAlignment','center','Rotation',90)
%         end
        
        if hh > 1
            title(titles{hh});
        end

        if ii == 3
            ylim([0 70]);
        end
        
        if hh == 1
            figPath = sprintf('V:\\retina\\John B\\phd backup\\papers\\Frontiers Degeneration of Visual Processing CFP\\img\\fig%d_%s',figNums(ii),figSuffixes{ii});
            export_fig(figPath,'-eps','-png','-transparent','-painters',pubfig);
%             print(gcf,figPath,'-depsc','-r600','-painters');
%             print(gcf,figPath,'-dpng','-r600','-painters');
            saveas(pubfig,figPath,'fig');

%             figPath = sprintf('V:\\retina\\John B\\phd backup\\presentations\\seattle jun 15\\fig%d_%s',figNums(ii),figSuffixes{ii});
%     
%             export_fig(figPath,'-m4','-png','-transparent','-painters');
        end
    end
    
    if hh == 1
        close(pubfig);
        figure(subfig);
    end     
end

%%

figure;
lineStyles = {'-' [] [] [] '--'};

for hh = 1:3
    subplot(1,3,hh);
    hold on;
    
    if hh == 1
        ylabel('Signal-to-Noise Ratio');
    end
    
    for ii = [1 5]
    %     ydata = cellfun(@(X) 100*X(:,:,ii)./repmat(X(:,1,1),1,6),snrData,'UniformOutput',false);

        for jj = 1:3
            [barhs(jj),m] = medianErrorbar([5 10 25 50 75 100],snrData{jj}(idxs{jj,hh},:,ii),'Color',colours(jj),'LineStyle',lineStyles{ii});
        end

        if ii == 1
            barhs(4) = line([-1 -2],[-1 -2],'Color','k','LineStyle','-');
            barhs(5) = line([-1 -2],[-1 -2],'Color','k','LineStyle','--');
            legend(barhs,{'18BGA' 'Flu' 'MFA' 'Control' 'Drug'});
        end
    end
    
    title(titles{hh});
end
    
xlabel('Flash Duration');
xlim([0 110]);

%%
    
medInfo = cell2mat(cellfun(@median,allInfoData,'UniformOutput',false));
uqInfo = cell2mat(cellfun(@(x) prctile(x,75),allInfoData,'UniformOutput',false))-medInfo;
lqInfo = medInfo-cell2mat(cellfun(@(x) prctile(x,25),allInfoData,'UniformOutput',false));

figure;

for ii = 1:2
    subplot(1,2,ii);
    barhs = barwitherr(cat(3,lqInfo(:,:,ii),uqInfo(:,:,ii)),medInfo(:,:,ii));
    barhs(end+1) = line(xlim,[12.5 12.5],'Color','k','LineStyle','--'); %#ok<SAGROW>
    legend(barhs,{'Control' 'Drug' 'Chance'});
    set(gca,'XTickLabel',uniqueDrugs);
    title(sprintf('Speed = %d{\\mu}m/s',1250/ii));

    if ii == 1
        ylabel('Mutual information (bits)');
    end
end

figFile = 'allbarinfo';
saveas(gcf,figFile,'fig');
saveas(gcf,figFile,'png');
close(gcf);

%%

baridxs = nsidxs;
baridxs(:,1) = cellfun(@(x) 1:size(x,1),perfData,'UniformOutput',false);
baridxs{1,2} = [1;3];
baridxs{1,3} = [2;4];
baridxs{3,2} = [1;2];
baridxs{3,3} = [3;4];
baridxs{1,4} = baridxs{1,1}(2:end);

%% MOVING BAR PERFORMANCE FIGURES
    
subfig = figure;
pubfig = figure;
set(pubfig,'Position',[100 100 780 400]);

for hh = 1:4
    medPerf = 100*cell2mat(cellfun(@(x,idx) median(x(idx,:,:)),perfData,baridxs(:,hh),'UniformOutput',false));
    uqPerf = 100*cell2mat(cellfun(@(x,idx) prctile(x(idx,:,:),75),perfData,baridxs(:,hh),'UniformOutput',false))-medPerf;
    lqPerf = medPerf-100*cell2mat(cellfun(@(x,idx) prctile(x(idx,:,:),25),perfData,baridxs(:,hh),'UniformOutput',false));
    
    for ii = 1:2
        if hh == 1
            subplot('Position',[0.075+0.5*(ii-1) 0.075 0.4 0.85]);
        else
            subplot(2,3,(ii-1)*3+hh-1);
        end
        
        [barhs,errhs] = barwitherr(cat(3,lqPerf(:,:,ii),uqPerf(:,:,ii)),medPerf(:,:,ii),'LineWidth',1.5);
        barhs(end+1) = line(xlim,[12.5 12.5],'Color','k','LineStyle','--','LineWidth',1.5); %#ok<SAGROW>
        
        for jj = 1:numel(errhs)
            set(errhs(jj),'LineWidth',1.5);
        end
        
        if ii == 2
            legend(barhs,{'Control' 'Drug' 'Chance'},'Location','North');
        end
        
        set(gca,'LineWidth',1.5,'XTickLabel',uniqueDrugs);
        title(sprintf('Speed = %d{\\mu}m/s - %s',1250/ii,titles{hh}));
        xlim([0.5 3.5])
        ylabel('Decoder Performance (%)');
        
        if hh == 1
            text(0,102.5,char('A'+ii-1),'FontSize',24,'FontWeight','bold');
        end
    end
    
    if hh == 1
        figure(subfig);
    end
end

%% SAVE MOVING BAR FIG FOR PRES

figPath = 'V:\\retina\\John B\\phd backup\\presentations\\seattle jun 15\\fig9_bar_perf';
    
export_fig(figPath,'-m4','-png','-transparent','-painters',pubfig);

%% SAVE MOVING BAR FIG FOR PUB

figPath = 'V:\\retina\\John B\\phd backup\\papers\\Frontiers Degeneration of Visual Processing CFP\\img\\fig9_bar_perf';
    
export_fig(figPath,'-eps','-png','-transparent','-painters',pubfig);
saveas(pubfig,figPath,'fig');
close(pubfig);