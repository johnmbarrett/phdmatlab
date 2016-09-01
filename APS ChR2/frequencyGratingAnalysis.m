load('channelNames.mat');
load('EventNo.mat');
load('spiketimestamps.mat');
load('contrast sequence.mat');

%%

ctrl = load('ChR2_cells_contrast_ctrl.mat','bestUnits');
drug = load('ChR2_cells_contrast_drug.mat','bestUnits');

commonUnits = intersect(ctrl.bestUnits,drug.bestUnits);
nCells = numel(commonUnits);

%%

ctrlRecording = 2;
drugRecording = 5;

gratingStimuli = find(~ismember(conditionOrder,[0 size(conditions,1)+1]));
gratingOnsets = 4*(gratingStimuli-1)+1;

stimulusTimes = [EventNo{ctrlRecording}(gratingOnsets); EventNo{drugRecording}(gratingOnsets)];
stimulusConditions = [repmat(conditions(conditionOrder(gratingStimuli),2:3),2,1) kron([0;1],ones(numel(gratingStimuli),1))];

%%

nSpikesss = cell(nCells,1);
nFigs = max(conditions(:,3))*2;
% figs = zeros(nFigs,1);

% for ii = 1:nFigs
%     figs(ii) = figure('Visible','off');
% end

for ii = 1:nCells
    tic;
    [~,~,~,~,~,~,nSpikesss{ii}] = rasterPlot(NaN,spiketimestamps{commonUnits(ii)},stimulusTimes,stimulusConditions(:,1),[-0.25 0 0.25 0.75],[],true,0.05,stimulusConditions(:,2:3),{'Phase' 'Frequency' 'Drug'});
    toc;
    continue;
    
    unitName = channelNames{1,commonUnits(ii)};
    
    for jj = 1:nFigs
        isCtrl = mod(jj,2) == 1;
        
        if isCtrl
            recording = ctrlRecording;
        else
            recording = drugRecording;
        end
        
        frequency = ceil(jj/2);
        
        figFile = sprintf('Stim %d\\freqraster_channel_%s_cluster_%s_freq_%d_drug_%d',recording,unitName(3:7),unitName(8),frequency,1-isCtrl);
        saveas(figs(jj),figFile,'fig');
        saveas(figs(jj),figFile,'png');
    end
    toc;
end

%%

temp = cell2mat(reshape(cellfun(@(A) sum(A(:,6:10),2),cat(3,nSpikesss{:}),'UniformOutput',false),[1 8 12 nCells]));
gratingResponses = zeros(50,8,6,2,nCells);

for ii = 1:2
    gratingResponses(:,:,:,ii,:) = temp(:,:,ii:2:end,:);
end

%%

% rotfun = @(pd) mod((0:7)+pd-4,8)+1;
% rotatedResponses = zeros(size(gratingResponses));
% 
% for ii = 1:nCells
%     for jj = 1:2
%         for kk = 1:6
%             [~,pd] = max(mean(gratingResponses(:,:,kk,jj,ii)));
%             rotatedResponses(:,rotfun(pd),kk,jj,ii) = gratingResponses(:,:,kk,jj,ii);
%         end
%     end
% end

%%

% maxResponses = mean(rotatedResponses(:,4,:,:,:));
% good = squeeze(all(all(maxResponses > 0)));
% normResponses = rotatedResponses(:,:,:,:,good)./repmat(maxResponses(:,:,:,:,good),[50 8 1 1 1]);
% meanResponses = squeeze(mean(normResponses));
% colours = 'rb';
% 
% for kk = 1:2
%     figure;
% 
%     for jj = 1:6
%         subplot(2,3,jj);
%         hold on;
% 
%     %     for kk = 1:2
%             errorbar(mean(meanResponses(:,jj,kk,:),4),std(meanResponses(:,jj,kk,:),[],4),'Color',colours(kk));
%     %     end
%     end
% end

%% 

mutualInformation = zeros(6,2,nCells);

x = kron((1:8)',ones(50,1));

for ii = 1:nCells
    tic;
    for jj = 1:2
        for kk = 1:6
            y = reshape(gratingResponses(:,:,kk,jj,ii),400,1);
            mutualInformation(kk,jj,ii) = biasCorrect(x,y);
        end
    end
    toc;
end

%%

Imin = lums([1 1 1:4]);
Imax = lums(5:end);
contrast = (Imax-Imin)./(Imax+Imin);

figure;
errorbar(repmat(barWidths',1,2),mean(mutualInformation,3),std(mutualInformation,[],3));
xlabel('Bar Width (pixels)');
ylabel('Mutual Information (nats)');
legend({'Control' 'MFA'});

%%

save('frequency_responses.mat','nSpikesss','gratingResponses','mutualInformation');
saveas(gcf,'info_vs_cont','fig');
saveas(gcf,'info_vs_cont','png');
close all;

%%

load('frequency_responses.mat');

%%

pid = pairwisePID(x,permute(reshape(gratingResponses,[400 6 2 nCells]),[4 1 2 3]));

%%

decoderPerformance = zeros(6,2);
perfVTrials = cell(6,2);
perfVUnits = cell(6,2);

% for ii = 1:nCells
%     tic;
    for jj = 1:2
        for kk = 1:6
            tic;
            y = reshape(gratingResponses(:,:,kk,jj,:),400,nCells);
            [fc,pvt,pvu] = simpleBayesianDecoder(x,y,5:5:40,2.^(0:7),10);
            decoderPerformance(kk,jj) = fc;
            perfVTrials{kk,jj} = pvt;
            perfVUnits{kk,jj} = pvu;
            toc;
        end
    end
%     toc;
% end

%%

figure;
errorbar(repmat(barWidths',1,2),mean(decoderPerformance,3),std(decoderPerformance,[],3));
xlabel('Bar Width (pixels)');
ylabel('Fraction Correct');
legend({'Control' 'MFA'});

%%

save('frequency_responses.mat','decoderPerformance','-append');
saveas(gcf,'freq_vs_perf','fig');
saveas(gcf,'freq_vs_perf','png');
close all;