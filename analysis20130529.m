dirs = {'F:\John B\Electrophys\JBWT0024', 'F:\John B\Electrophys\JBWT0025', 'D:\temp\project\retina data\Wild-Type Stim\JBWT0026', 'C:\temp\retina data\JBWT0027', 'D:\temp\project\retina data\Wild-Type Stim\JBWT0028', 'D:\temp\project\retina data\Wild-Type Stim\JBWT0029'};
sigs = {'rf', 'visual', 'visual' 'visual', 'visual', 'visual'};
recs = {[2 4 7 9], [2 4], 4, [7 9], [2 4], [2 4]};

assert(numel(dirs) == numel(sigs) && numel(sigs) == numel(recs));
currentDir = pwd;

titles = {'P26 Retina A' 'P26 Retina B' 'P99 Retina A' 'P12 Retina A' 'P13 Retina B' 'P18 Retina A' 'P73 Retina A'};

%%
% for ii = 1:numel(dirs)
%     rec = recs{ii};
%     
%     cd(dirs{ii});
%     recordings = initRecordings;
%     
%     for jj = 1:numel(rec)
%         try
%             recording = recordings(rec(jj));
%             analyseGollischMeisterResponses(recording,'latencytype','peak','significance',sigs{ii});
%             analyseJacobsNirenbergResponses(recording,'latencytype','peak','significance',sigs{ii});
%         catch err
%             logMatlabError(err);
%             cd(currentDir);
%         end
%     end
% end
% 
% cd(currentDir);

%%

figure

subplots = 1;
colours = [1 0 0; 0 1 1; 0 0 1; 1 0 1];
ylims = [0 0.6; 0 0.7; -0.1 0.8; 0 0.9; 0 0.9; 0 0.4; -0.1 0.7];

for ii = 1:numel(dirs);
    cd(dirs{ii});
    
    try
        recordings = initRecordings;

        subplot(4,2,subplots);
        hold on;

        for jj = 1:numel(recs{ii});
            if jj == 3
                subplots = subplots + 1;
                subplot(4,2,subplots);
                hold on;
            end

            rec = recs{ii}(jj);
            recording = recordings(rec);
            outDir = getAnalysisOutputDir(recording);
            load([outDir '\' recording.dataFile '_info']);
            load([outDir '\' recording.dataFile '_psrh_abs'],'valuess');
            count = mutualCountPhaseInformationPerWidth2(:,:,:,2);
            lat = mutualLatencyPhaseInformationPerWidth4;
            isLight = 1-mod(floor(rec/2),2);
            x = sort(25*valuess{1}/2);
            errorbar(x,mean(count(:,:,2)),2*std(count(:,:,2))/sqrt(size(count,1)),'Color',colours(3-isLight,:),'LineStyle','--');
            errorbar(x,mean(lat(:,:,2)),2*std(lat(:,:,2))/sqrt(size(lat,1)),'Color',colours(3*isLight+1,:),'LineStyle','--');
            errorbar(x,mean(count(:,:,1)),2*std(count(:,:,1))/sqrt(size(count,1)),'Color',colours(3-isLight,:));
            errorbar(x,mean(lat(:,:,1)),2*std(lat(:,:,1))/sqrt(size(lat,1)),'Color',colours(3*isLight+1,:));
            
            if isLight
                set(gca,'XTick',x);
                title(titles{subplots});
                xlabel('Bar Width/{\mu}m');
                xlim([x(1)-50 x(end)+50]);
                ylabel('Mutual information/bits');
                ylim(ylims(subplots,:));
%                 ylim([0 1.1*max(max(mean(count)+std(count)/sqrt(size(count,1))))]);
            end
        end

        subplots = subplots + 1;
    catch err
        logMatlabError(err);
        cd(currentDir);
        return;
    end
end
    
cd(currentDir);

%%

figure

subplots = 1;
colours = 'rgbmyc';
ylims = [0 0.7; 0 0.7; 0 0.45; 0 0.45; 0 0.5; 0 0.62; 0 0.7];

for ii = 1:numel(dirs);
    cd(dirs{ii});
    
    try
        recordings = initRecordings;

        subplot(4,2,subplots);
        hold on;

        for jj = 1:numel(recs{ii});
            if jj == 3
                subplots = subplots + 1;
                subplot(4,2,subplots);
                hold on;
            end

            rec = recs{ii}(jj);
            recording = recordings(rec);
            outDir = getAnalysisOutputDir(recording);
            load([outDir '\' recording.dataFile '_bayes']);
            load([outDir '\' recording.dataFile '_psrh_abs'],'valuess');
            count = countPerformance;
            time = timingPerformance;
            corr = correlationPerformance;
            isLight = 1-mod(floor(rec/2),2);
            x = sort(25*valuess{1}/2);
            errorbar(x,mean(count(:,:,2)),2*std(count(:,:,2))/sqrt(size(count,1)),'Color',colours(3*isLight+3),'LineStyle','--');
            errorbar(x,mean(count(:,:,1)),2*std(count(:,:,1))/sqrt(size(count,1)),'Color',colours(3*isLight+3));
            errorbar(x,mean(time(:,:,2)),2*std(time(:,:,2))/sqrt(size(time,1)),'Color',colours(3*isLight+1),'LineStyle','--');
            errorbar(x,mean(time(:,:,1)),2*std(time(:,:,1))/sqrt(size(time,1)),'Color',colours(3*isLight+1));
            errorbar(x,mean(corr(:,:,2)),2*std(corr(:,:,2))/sqrt(size(corr,1)),'Color',colours(3*isLight+2),'LineStyle','--');
            errorbar(x,mean(corr(:,:,1)),2*std(corr(:,:,1))/sqrt(size(corr,1)),'Color',colours(3*isLight+2));
            line([x(1)-50 x(end)+50],[1 1]/8,'Color','k','LineStyle',':');
            
            if isLight
                set(gca,'XTick',x);
                title(titles{subplots});
                xlabel('Bar Width/{\mu}m');
                xlim([x(1)-50 x(end)+50]);
                ylabel('Fraction Correct');
                ylim(ylims(subplots,:));
%                 ylim([0 1.1*max(max(mean(count)+std(count)/sqrt(size(count,1))))]);
            end
        end

        subplots = subplots + 1;
    catch err
        logMatlabError(err);
        cd(currentDir);
        return;
    end
end
    
cd(currentDir);