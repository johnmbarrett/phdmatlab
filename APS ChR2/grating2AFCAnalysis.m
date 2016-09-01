function grating2AFCAnalysis(independentVariable,ctrlRecording,drugRecording,useAllResponders,blocks,decoder,parallelisation,ctrlCellsSuffix,drugCellsSuffix,saveFileSuffix)
% confounder = 'contrast';
% confounder = 'frequency';
    if nargin < 10
        saveFileSuffix = '';
    else
        saveFileSuffix = ['_' saveFileSuffix];
    end
    
    if nargin < 9
        drugCellsSuffix = 'drug';
    end
    
    if nargin < 8
        ctrlCellsSuffix = 'ctrl';
    end

    if nargin < 7
        parallelisation = struct('section','none');
    end
    
    if nargin < 6
        decoder = 'discrete';
    end

    if nargin < 4
        useAllResponders = false;
    end

    %%
    
    load(['./' independentVariable ' sequence.mat']);

    %%

    if nargin < 5
        blocks = (1:nReps)';
    elseif isscalar(blocks)
        blocks = (1:blocks)';
    else
        blocks = blocks(:);
    end

    nBlocks = numel(blocks);
    
    %%

    if strcmp(independentVariable,'contrast')
        legendLocation = 'NorthWest';
        xvals = 100*C';
        xticks = xvals;
        xticklabels = arrayfun(@(c) sprintf('%d',c),10*floor(xvals/10),'UniformOutput',false);
        xlab = 'Michelson Contrast (%)';
        xlims = [0 70];
        xscale = 'linear';
%         splits = [1 3; 4 6];
%         subplotTitles = arrayfun(@(c) sprintf('Contrast = %3.1f%%',c),100*C,'UniformOutput',false);
    elseif strcmp(independentVariable,'frequency')
        legendLocation = 'NorthEast';
        xvals = (30./(barWidths*4*2))';
        xticks = wrev(xvals);
        xticklabels = arrayfun(@(f) sprintf('%0.3f',f),xticks,'UniformOutput',false);
        xlab = 'Frequency (cpd)';
        xlims = [10^-2 10^-0.4];
        xscale = 'log';
%         splits = [1 6];
%         subplotTitles = arrayfun(@(f) sprintf('Frequency = %0.4f cpd',f),xvals,'UniformOutput',false);
    else
        legendLocation = 'NorthWest';
        xvals = 100*(Lmed-Lmin)'./Lmin;
        xticks = xvals;
        xticklabels = arrayfun(@(l) sprintf('%4.0f',l),xvals,'UniformOutput',false);
        xlab = 'Weber Contrast (%)';
        xlims = [0 1500];
        xscale = 'linear';
    end
    
    doRastersOnly = any(strcmp(parallelisation.section,'rasters'));
    doExcitingPlots = any(strcmp(parallelisation.section,'plots'));
    doInit = doRastersOnly || doExcitingPlots || any(strcmp(parallelisation.section,{'none' 'init'}));
    doMain = ~doRastersOnly && doExcitingPlots && any(strcmp(parallelisation.section,{'none' 'main'}));
    doFinish = ~doRastersOnly && doExcitingPlots && any(strcmp(parallelisation.section,{'none' 'finish'}));
    
    initFile = sprintf('gratings_perf_vs_%s_decoder_%s_init%s.mat',independentVariable,decoder,saveFileSuffix);
    likelihoodsFilePrefix = sprintf('gratings_perf_vs_%s_decoder_%s_likelihoods%s',independentVariable,decoder,saveFileSuffix);
    
    resultsFile = [independentVariable '_2afc_responses_decoder_' decoder saveFileSuffix '.mat'];
    
    %%

    if doInit
        %%
        load('channelNames.mat');
        load('EventNo.mat');
        load('spiketimestamps.mat');

        %%

        ctrl = load(['ChR2_cells_' independentVariable '_' ctrlCellsSuffix '.mat'],'bestUnits');
        drug = load(['ChR2_cells_' independentVariable '_' drugCellsSuffix '.mat'],'bestUnits');
        units = cell(1,2);

        if ~islogical(useAllResponders)
            niceUnits = {39 42 'a'};
    %             5	29	'a';	...
    %                      6	23	'a';	...
    %                      9	19	'a';	...
    %                      12	16	'a';	...
    %                      23	28	'a';	...
    %                      26	22	'a';	...
    %                      28	32	'b';	...
    %                      29	30	'b';	...
    %                      29	32	'b';	...
    %                      30	11	'a';	...
    %                      30	30	'a';	...
    %                      31	45	'a';	...
    %                      32	31	'a';	...
    %                      33	48	'a';	...
    %                      33	57	'a';	...
    %                      38	43	'b';	...
    %                      39	42	'a'};

             niceUnitNames = cellfun(@(r,c,u) sprintf('Ch%02d_%02d%c',r,c,u),niceUnits(:,1),niceUnits(:,2),niceUnits(:,3),'UniformOutput',false);

             units{1} = find(ismember(channelNames(1,:)',niceUnitNames)); %#ok<NODEF>
             units{2} = units{1};
        elseif useAllResponders
            allUnits = union(ctrl.bestUnits,drug.bestUnits);

            for stimulusIndex = 1:2
                units{stimulusIndex} = allUnits;
            end
        else
            units{1} = ctrl.bestUnits;
            units{2} = drug.bestUnits;
        end
        
        nCells = cellfun(@numel,units);

        %%

    %     ctrlRecording = 2;
    %     drugRecording = 5;
        recordings = [ctrlRecording drugRecording];

        gratingStimuli = find(~ismember(conditionOrder,[0 size(conditions,1)+1]));
        nConditions = size(conditions,1);
        gratingStimuli = gratingStimuli(kron(blocks-1,ones(nConditions,1))*nConditions+repmat((1:nConditions)',nBlocks,1));
        stimulusConditions = conditions(conditionOrder(gratingStimuli),2:4);
        trialOrder = stimulusConditions(:,3);

        firstStimOnsets = 6*(gratingStimuli-1)+1;
    %     secondStimOnsets = 6*(gratingStimuli-1)+2;
        thirdStimOnsets = 6*(gratingStimuli-1)+3;

        gratingOnsets = zeros(size(firstStimOnsets));
        gratingOnsets(trialOrder == 0) = firstStimOnsets(trialOrder == 0);
        gratingOnsets(trialOrder == 1) = thirdStimOnsets(trialOrder == 1);

        maskOnsets = zeros(size(firstStimOnsets));
        maskOnsets(trialOrder == 0) = thirdStimOnsets(trialOrder == 0);
        maskOnsets(trialOrder == 1) = firstStimOnsets(trialOrder == 1);

        gratingTimes = cell2mat(cellfun(@(t) t(gratingOnsets),EventNo(recordings),'UniformOutput',false)); %#ok<USENS>
        maskTimes = cell2mat(cellfun(@(t) t(maskOnsets),EventNo(recordings),'UniformOutput',false));

        s = stimulusConditions(:,1)+4*(stimulusConditions(:,2)-1);

        %%

        [uniqueConditions,~,phaseXLevel] = unique(stimulusConditions(:,1:2),'rows');
        disp(sum(cellfun(@numel,units)));
        
        %%

        if doRastersOnly
            for ii = 1:2
                for jj = 1:numel(units{ii})
                    tic;
                    figs = rasterPlot(spiketimestamps{units{ii}(jj)},EventNo{recordings(ii)}(firstStimOnsets),phaseXLevel,-0.25:0.25:1,[],true,0.05,stimulusConditions(:,3),{'Phase X Level','Mask First'},true); %#ok<USENS>

                    for kk = 1:2
                        figFile = sprintf('Stim %d\\2afcraster_channel_%s_maskfirst_%d',recordings(ii),channelNames{1,units{ii}(jj)},kk-1); %#ok<USENS>
                        saveas(figs(kk),figFile,'fig');
                        saveas(figs(kk),figFile,'png');
                        close(figs(kk));
                    end

                    toc;
                end
            end
            
            return;
        end
        
        if doExcitingPlots
            binWidth = 0.001;
%             bandWidth = 0.01;
%             kernel = normpdf(-5*bandWidth/binWidth:5*bandWidth/binWidth,0,bandWidth/binWidth);
            growthConstant = 0.001;
            decayConstant = 0.02;
            t = 0:5*(max(growthConstant,decayConstant)/binWidth);
            psp = (1-exp(-t/(growthConstant/binWidth))).*exp(-t/(decayConstant/binWidth));
            kernel = [zeros(1,t(end)) psp];
            kernel = kernel/sum(kernel);
                
            timeIndices = -0.25:binWidth:1;
            timeIndices = timeIndices(ceil(numel(kernel)/2):end-floor(numel(kernel)/2));
            
            normDGratingSpikess = cell(1,2);
            normDGratingFRs = cell(1,2);
            
            for ii = 1:2
                nCells = numel(units{ii});
                squareHist = zeros(nCells,nReps,1.25/binWidth+1,4,6,2);
                
                for jj = 1:nCells
                    tic;
                    [~,lines] = rasterPlot(NaN,spiketimestamps{units{ii}(jj)},EventNo{recordings(ii)}(firstStimOnsets),phaseXLevel,-0.25:0.25:1,[],true,0.05,stimulusConditions(:,3),{'Phase X Level','Mask First'},true); %#ok<USENS>
                    
                    for kk = 1:24
                        phase = uniqueConditions(kk,1);
                        level = uniqueConditions(kk,2);
                        
                        for ll = 1:nReps
                            for mm = 1:2
                                binIndices = ceil((lines{kk,mm}{ll}+0.25)/binWidth);
                                squareHist(jj,ll,binIndices,phase,level,mm) = 1;
                            end
                        end
                    end
                    
                    for kk = 1:6
                        gratingFirstResponse = sum(sum(squareHist(jj,:,(0.25/binWidth+1):(0.5/binWidth),:,kk,1),2),3);
                        gratingLastResponse = sum(sum(squareHist(jj,:,(0.75/binWidth+1):(1/binWidth),:,kk,1),2),3);
                        gratingResponse = squeeze(gratingFirstResponse + gratingLastResponse);
                        [~,bestPhase] = max(gratingResponse);
                        squareHist(jj,:,:,:,kk,:) = squareHist(jj,:,:,circshift(1:4,[0 2-bestPhase]),kk,:);
                    end
                    toc;
                end
                
                squareHist2 = [squareHist(:,:,:,:,:,1) squareHist(:,:,[ceil(size(squareHist,3)/2):size(squareHist,3) 1:floor(size(squareHist,3)/2)],:,:,2)];
                dGratingSpikes = squeeze(diff(mean(sum(squareHist2(:,:,(0.25/binWidth+1):0.5/binWidth,[4 2],:),3),2),[],4));
                maskSpikes = mean(reshape(sum(squareHist2(:,:,(0.75/binWidth+1):1/binWidth,:,:),3),[nCells nReps*4*6*2]),2);
                normDGratingSpikess{ii} = dGratingSpikes./repmat(maskSpikes,[1 6]);

                squareHist3 = squeeze(mean(squareHist,1));
                frs = zeros(nReps,size(squareHist3,2)-numel(kernel)+1,4,6,2);
                
                for jj = 1:2
                    for kk = 1:6
                        for ll = 1:4
                            for mm = 1:nReps
                                frs(mm,:,ll,kk,jj) = conv(squareHist3(mm,:,ll,kk,jj),kernel,'valid');
                            end
                        end
                    end
                end
                
                frs = [frs(:,:,:,:,1); frs(:,[ceil(size(frs,2)/2):size(frs,2) 1:floor(size(frs,2)/2)],:,:,2)];
                
                dGratingFR = squeeze(diff(trapz(binWidth:binWidth:0.25,frs(:,floor(numel(kernel)/2)+(1:0.25/binWidth),[4 2],:),2),[],3));
                maskFR = mean(reshape(trapz(binWidth:binWidth:0.25,frs(:,floor(numel(kernel)/2)+0.5/binWidth+(1:0.25/binWidth),:,:),2),nReps*2*4*6,1));
                normDGratingFRs{ii} = dGratingFR./repmat(maskFR,50,6);
                
                figure;
                set(gcf,'Position',[100 100 1000 750]);
                axs = zeros(6,1);
                ylimits = [Inf -Inf];
                
                for jj = 1:6
                    axs(jj) = subplot(2,3,jj);
                    boundedline(timeIndices,squeeze(mean(frs(:,:,:,jj)/binWidth)),permute(std(frs(:,:,:,jj)/binWidth),[2 1 3]));
                    yy = ylim;
                    ylimits(1) = min(ylimits(1),yy(1));
                    ylimits(2) = max(ylimits(2),yy(2));
                end
                
                for jj = 1:6
                    hold(axs(jj),'on');
                    
                    plot(axs(jj),repmat(0:0.25:0.75,2,1),repmat(ylimits',1,4),'Color','k','LineStyle','--');
                    
                    set(axs(jj),'Position',[0.07+0.32*mod(jj-1,3) 0.54-0.48*(jj>3) 0.27 0.41]);
                    
                    title(axs(jj),sprintf('Frequency = %0.3f cpd',xvals(jj)))
                    xlim(axs(jj),timeIndices([1 end]));
                    ylim(axs(jj),ylimits);
                    
                    if jj == 4
                        xlabel(axs(jj),'Time (s)');
                        ylabel(axs(jj),'Firing Rate (Hz)');
                    end
                end
                
                figFile = sprintf('2afcgratings_sdf_drug_%d',ii);
                saveas(gcf,figFile,'fig');
                export_fig(figFile,'-eps','-png','-m4','-transparent','-painters');
                close(gcf);
            end
            
            figure;
            hold on;
            colourOrder = get(gca,'ColorOrder');
            
            for ii = 1:2
                p = squeeze(prctile(normDGratingSpikess{ii},[25 50 75],1));
                m = p(2,:)';
                l = m-p(1,:)';
                u = p(3,:)'-m;
                errorbar(xvals,m,l,u,'Color',colourOrder(ii,:));
            end
            
            set(gca,'XScale','log');
            legend({'Control' 'MFA'});
            xlabel('Spatial Frequency (cpd)');
            xlim(10.^[-2 -0.4])
            ylabel('(Preferred Grating Spikes-Null Grating Spikes)/Mask Spikes');
            
            figFile = '2afcgratings_dspikes_vs_freq';
            saveas(gcf,figFile,'fig');
            export_fig(figFile,'-eps','-png','-m4','-transparent','-painters');
            close(gcf);
            
            figure;
            hold on;
            
            for ii = 1:2
                p = prctile(normDGratingFRs{ii},[25 50 75]);
                m = p(2,:)';
                l = m-p(1,:)';
                u = p(3,:)'-m;
                errorbar(xvals,m,l,u,'Color',colourOrder(ii,:));
            end
            
            set(gca,'XScale','log');
            legend({'Control' 'MFA'});
            xlabel('Spatial Frequency (cpd)');
            xlim(10.^[-2 -0.4])
            ylabel('(Preferred Grating FR-Null Grating FR)/Mask FR');
            
            figFile = '2afcgratings_dfr_vs_freq';
            saveas(gcf,figFile,'fig');
            export_fig(figFile,'-eps','-png','-m4','-transparent','-painters');
            close(gcf);
            
            if exist(resultsFile,'file')
                save(resultsFile,'-append','normDGratingSpikess','normDGratingFRs');
            else
                save(resultsFile,'normDGratingSpikess','normDGratingFRs');
            end
            
            return;
        end

        %%

        Gs = cell(1,2);
        Ms = cell(1,2);

        %%

        % this section of code makes Baby Knuth cry :'(
    %     units2 = {2559 2559};
        if strcmp(decoder,'discrete')
            uniformOutput = true;
            respFun = @(st,t0) sum(st > t0 & st <= t0 + 0.25);
            catFun = @cell2mat;
            initFun = @initDiscreteDecoder;
            trainFun = @trainDiscreteDecoder;
            testFun = @getDiscreteLikelihood;
        elseif strncmp(decoder,'jacobs',6)
            uniformOutput = false;
    %         respFun = @(st,t0) st(st > t0 - 0.25 & st <= t0 + 0.5) - t0;
            respFun = @(st,t0) st(st > t0 & st <= t0 + 0.25) - t0;
            catFun = @(Cs) [Cs{:}];

            if numel(decoder) > 6
                densityEstimationMethod = decoder(7:end);
            else
                densityEstimationMethod = 'kde';
            end

            [trainFun,initFun] = getJacobsDecoderTrainingFun(0.001,0.005,0.25,densityEstimationMethod);
            testFun = @getJacobsLikelihood;
        else
            error('Unknown decoder');
        end

        %%

        for hh = 1:2
            Gs{hh} = catFun(cellfun(@(st) arrayfun(@(gt) respFun(st,gt),gratingTimes(:,hh),'UniformOutput',uniformOutput),spiketimestamps(units{hh}),'UniformOutput',false));
            Ms{hh} = catFun(cellfun(@(st) arrayfun(@(mt) respFun(st,mt),maskTimes(:,hh),'UniformOutput',uniformOutput),spiketimestamps(units{hh}),'UniformOutput',false));
        end
        
        %%
        
%         normResp = cellfun(@(u) zeros(nReps,numel(u),2,4,6,2),units,'UniformOutput',false);
%         maxResp = cellfun(@(G,M) repmat(max([G;M],[],1),nReps,1),Gs,Ms,'UniformOutput',false);
%         
% %         assert(all(cellfun(@(A) all(A(:) > 0),maxResp)));
%         
%         %%
%         
%         for ii = 1:2
%             for jj = 1:6
%                 for kk = 1:4
%                     for ll = 1:2
%                         stimIndices = stimulusConditions(:,1) == kk & stimulusConditions(:,2) == jj & stimulusConditions(:,3) == ll-1;
%                         normResp{ii}(:,:,ll,kk,jj,1) = Gs{ii}(stimIndices,:)./maxResp{ii};
%                         normResp{ii}(:,:,ll,kk,jj,2) = Ms{ii}(stimIndices,:)./maxResp{ii};
%                     end
%                 end
%                 
%                 for kk = 1:nCells(ii)
%                     phaseResp = permute(mean(mean(mean(normResp{ii}(:,kk,:,:,:,1),1),3),5),[4 1 2 3 5 6]);
%                     bestPhase = find(phaseResp == max(phaseResp,1));
%                     normResp{ii} = circshift(normResp{ii},[0 0 0 2-bestPhase 0 0]);
%                 end
%             end
%         end
%         
%         %%
%         
%         dResp = reshape(mean(normResp{2}(:,:,:,:,:,1),1),nCells(2),48);
%         dResp = max(dResp,[],2)-min(dResp,[],2);
%         [dResp,sortIndices] = sort(dResp,'descend');
%         normResp{2} = normResp{2}(:,sortIndices,:,:,:,:);
%         
%         %%
%         
%         figure;
%         set(gcf,'Position',[9 49 944 948]);
%         
%         for ii = 1:6
%             subplot('Position',[0.06+0.49*mod(ii-1,2) 1.0225-0.31*ceil(ii/2) 0.44 0.25]);
%             errorbar(permute(mean(reshape(permute(normResp{2}(:,find(~isnan(dResp),1),:,:,ii,:),[1 3 4 6 2 5]),[50 4 2]),1),[2 3 1]),permute(std(reshape(permute(normResp{2}(:,find(~isnan(dResp),1),:,:,ii,:),[1 3 4 6 2 5]),[50 4 2]),[],1),[2 3 1]),'LineWidth',1.5);
% %             errorbar(reshape(permute(mean(permute(normResp{2}(:,find(~isnan(dResp),1),:,:,ii,:),[1 3 4 6 2 5]),1),[3 2 4 1 5 6]),[4 4 1 1 1 1]),reshape(permute(std(permute(normResp{2}(:,find(~isnan(dResp),1),:,:,ii,:),[1 3 4 6 2 5]),1),[3 2 4 1 5 6]),[4 4 1 1 1 1]),'LineWidth',1.5);
%             
%             set(gca,'LineWidth',1.5,'XTick',1:4);
%             
%             title(sprintf('Contrast = %d0%%',ii));
%             
%             if ii == 5
%                 xlabel('Phase');
%                 ylabel('Normalised Response');
%             end
%             
%             ylim([-0.2 0.8]);
%         end
%         
%         lh = legend({'Grating' 'Mask'},'LineWidth',1.5,'Orientation','Horizontal');
% %         lh = legend({'Grating (Grating First)' 'Grating (Mask First)' 'Mask (Grating First)' 'Mask (Mask First)'},'LineWidth',1.5,'Orientation','Horizontal');
%         set(lh,'Position',[0.06 0.01 0.93 0.03])
%         
%         %%
%         
%         figure;
%         set(gcf,'Position',[9 49 944 948]);
%         offset = 2;
%         
%         for ii = 1:2
%             subplot('Position',[0.06+0.49*mod(ii-1,2) 1-0.49*ceil(ii/2) 0.44 0.45]);
%             
% %             p = prctile(reshape(permute(normResp{2}(:,find(~isnan(dResp),1)+offset,:,ii,:,:),[1 3 5 6 2 4]),[50 6 2]),[25 50 75]);
% %             m = squeeze(p(2,:,:));
% %             l = m-squeeze(p(1,:,:));
% %             u = squeeze(p(3,:,:))-m;
% %             
% %             errorbar(repmat((0.1:0.1:0.6)',1,2),m,l,u,'LineWidth',1.5);
%             
%             errorbar(repmat((0.1:0.1:0.6)',1,2),permute(mean(reshape(permute(normResp{2}(:,find(~isnan(dResp),1)+offset,:,ii,:,:),[1 3 5 6 2 4]),[50 6 2]),1),[2 3 1]),permute(std(reshape(permute(normResp{2}(:,find(~isnan(dResp),1)+offset,:,ii,:,:),[1 3 5 6 2 4]),[50 6 2]),[],1),[2 3 1]),'LineWidth',1.5);
%             
%             set(gca,'LineWidth',1.5,'XTick',0:0.1:0.7,'XTickLabel',0:10:70);
%             
%             title(sprintf('Phase %d',ii));
%             
%             xlabel('Michelson Contrast (%)');
%             
%             if ii == 1
%                 legend({'Grating' 'Mask'},'LineWidth',1.5,'location','SouthWest');
%                 ylabel('Normalised Response');
%             end
%             
%             ylim([-0.05 0.55]);
%         end
%         
%         %%
%         
%         idx = find(~isnan(dResp),1)+offset;
%         idx = sortIndices(idx);
%         idx = units{2}(idx);
%         chlabel = channelNames{1,idx};
%         [x,y] = ellipseFn(str2double(chlabel(6:7))*42/4-21/4,64*42/4-str2double(chlabel(3:4))*42/4+21/4,250/4,250/4,0);
%         
%         %%
%         
% %         Is = {I1 I2};d
%         for ii = 3:4
%             subplot('Position',[0.06+0.49*mod(ii-1,2) 1-0.49*ceil(ii/2) 0.44 0.42]);
%             hold on;
% %             imshow(Is{ii-2});
% %             pcolor(flipud(double(Is{ii-2})));
% %             colormap(gray);
% %             caxis([0 255]);
% %             view(2);
% %             shading flat;
%             fill([0 664 664 0],[0 0 664 664],56*[1 1 1]/255);
%             
%             for jj = 1:4
%                 fill([0 664 664 0],664-[0 0 160 160]-320*(jj-1)-80*(ii-3)+1,248*[1 1 1]/255,'EdgeColor','none');
%             end
%             
%             ez = ezplot(x,y);
%             box on;
%             set(ez,'LineWidth',1.5);
%             
%             set(gca,'XTick',[],'YTick',[]);
%             title('');
%             xlabel('');
% %             xlim([0 664])
%             xlim([-570 570]/4+str2double(chlabel(6:7))*42/4-21/4);
%             
%             if ii == 3
%                 ylabel('Receptive Field');
%             else
%                 ylabel('');
%             end
%             
% %             ylim([0 664]);
%             ylim([-500 500]/4+64*42/4-str2double(chlabel(3:4))*42/4+21/4);
%             
%             fill(570*[-1 1 1 -1]/4+str2double(chlabel(6:7))*42/4-21/4+[1 0 0 1],500*[-1 -1 1 1]/4+64*42/4-str2double(chlabel(3:4))*42/4+21/4-[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'FaceColor','none','LineWidth',1.5);
%             set(gca,'Position',[0.035+0.49*mod(ii-1,2) 1-0.495*ceil(ii/2) 0.49 0.44],'XTick',[],'YTick',[]);  
%         end
%         
%         %%
%         
%         saveas(gcf,'V:\retina\John B\phd backup\thesis\img\apschr2\new_contrast_model_example','fig');
%         export_fig('V:\retina\John B\phd backup\thesis\img\apschr2\new_contrast_model_example','-eps','-png',-'m4','-transparent','-painters');
        
        %%
        
%         figure; hold on; imshow(I1); ezplot(x,y);
%         figure; hold on; imshow(I2); ezplot(x,y);
%         figure; hold on; imshow(I3); ezplot(x,y);
%         figure; hold on; imshow(I4); ezplot(x,y);
        
        %%
        
%         exampleGratResp = reshape(permute(normResp{2}(:,find(~isnan(dResp),1),:,:,:,1),[1 3:5 2 6]),[50 4 6])*maxResp{2}(1,sortIndices(find(~isnan(dResp),1))); %#ok<NASGU>
%         exampleMaskResp = reshape(permute(normResp{2}(:,find(~isnan(dResp),1),:,:,:,2),[1 3:5 2 6]),[1200 1])*maxResp{2}(1,sortIndices(find(~isnan(dResp),1))); %#ok<NASGU>
%         save('V:\retina\John B\phd backup\matlab\APS ChR2\contrast_example_data.mat','exampleGratResp','exampleMaskResp');

        %%

        nStimuli = size(gratingStimuli,1);
        t = kron([2;1],ones(nStimuli,1));
    %     t = [s+1; ones(nStimuli,1)];
    %     logPs = log([0.5 repmat(1/48,1,24)]);

        Rs = cellfun(@vertcat,Gs,Ms,'UniformOutput',false);
        paramss = cellfun(@(R) initFun(t,R),Rs,'UniformOutput',false);
        trials = repmat((1:nStimuli)',2,1)-1;
        
        save(initFile,'trials','nStimuli','t','Gs','Ms','Rs','paramss','trainFun','testFun','stimulusConditions');
    end
    
    %%
    
    if doMain
        load(initFile);
        
        if strcmp(parallelisation.section,'none')
            stimulusIter = 1:nStimuli;
            conditionIter = 1:2;
            cellIter = cellfun(@(R) 1:size(R,2),Rs,'UniformOutput',false);
        else
            stimulusIter = parallelisation.stimulusIter;
            conditionIter = parallelisation.conditionIter;
            cellIter = parallelisation.cellIter;
        end
            
%         detectionSuccess = zeros(nBlocks*2,4,6,2);
%         seen = zeros(4,6);

        gratingLikelihoods = cellfun(@(cells) zeros(numel(stimulusIter),numel(cells)),cellIter,'UniformOutput',false);
        maskLikelihoods = cellfun(@(cells) zeros(numel(stimulusIter),numel(cells)),cellIter,'UniformOutput',false);

    %     discriminationSuccess = zeros(nReps*2,4,6,2);

        for ii = 1:numel(stimulusIter)
            stimulusIndex = stimulusIter(ii);
            tic;
            testTrial = mod(trials,nStimuli) == stimulusIndex-1;
            gratingTest = find(testTrial(1:nStimuli));

            phase = stimulusConditions(gratingTest,1);
            level = stimulusConditions(gratingTest,2);
%             seen(phase,level) = seen(phase,level) + 1;

            gratingTrainIndices = find(~testTrial(1:nStimuli) & stimulusConditions(:,1) == phase & stimulusConditions(:,2) == level);

            if strcmp(independentVariable,'luminance')
                maskTrainIndices = find(~testTrial(nStimuli+1:2*nStimuli) & stimulusConditions(:,2) == level)+nStimuli;
            else
                maskTrainIndices = find(~testTrial(nStimuli+1:2*nStimuli))+nStimuli;
            end

            trainIndices = [gratingTrainIndices; maskTrainIndices];

            PDFs = cellfun(@(R,params,cells) trainFun(trainIndices,t,R(:,cells),params),Rs(conditionIter),paramss(conditionIter),cellIter);

            for jj = 1:numel(conditionIter)
                conditionIndex = conditionIter(jj);
                gratingResponse = Gs{conditionIndex}(gratingTest,cellIter{jj});
                gratingLikelihood = testFun(gratingResponse,PDFs(conditionIndex));
                gratingLikelihood(isinf(gratingLikelihood)) = log(eps);
%                 gratingLikelihood = sum(gratingLikelihood,1);
                gratingLikelihoods{jj}(ii,:) = gratingLikelihood(:,2);
%                 gratingPosterior = gratingLikelihood+logPs;

                maskResponse = Ms{conditionIndex}(testTrial(nStimuli+1:2*nStimuli),cellIter{jj});
                maskLikelihood = testFun(maskResponse,PDFs(conditionIndex));
                maskLikelihood(isinf(maskLikelihood)) = log(eps);
%                 maskLikelihood = sum(maskLikelihood,1);
                maskLikelihoods{jj}(ii,:) = maskLikelihood(:,2);
%                 maskPosterior = maskLikelihood+logPs;
            end
            toc;
        end
        
        if strcmp(parallelisation.section,'none')
            likelihoodsFile = likelihoodsFilePrefix;
        else
            likelihoodsFile = sprintf('%s_trials_%s',likelihoodsFilePrefix,arrayToFilenameString(stimulusIter));
        
            conditionNames = {'control' 'drug'};

            for ii = 1:2
                if ismember(ii,conditionIter)
                    likelihoodsFile = sprintf('%s_%s_cells_%s',likelihoodsFile,conditionNames{ii},arrayToFilenameString(cellIter{conditionIter == ii}));
                end
            end
        end
        
        save(likelihoodsFile,'gratingLikelihoods','maskLikelihoods','conditionIter','cellIter','stimulusIter');
    end
    
    %%

    if ~doFinish
        return;
    end
    
    if ~doInit
        load(initFile,'stimulusConditions');
    end
    
    load(sprintf('%s.mat',likelihoodsFilePrefix));
    
    detectionSuccess = zeros(nBlocks*2,4,6,2);
    logPs = [0.5 0.5];
    
    for ii = 1:2
        gratingPosterior = sum(gratingLikelihoods{ii},2)+logPs(2);
        maskPosterior = sum(maskLikelihoods{ii},2)+logPs(2);
        
        success = zeros(size(gratingPosterior));
        
        success(gratingPosterior > maskPosterior) = 1;
        success(gratingPosterior == maskPosterior) = 0.5;
    
        for jj = 1:6
            for kk = 1:4
                idx = stimulusConditions(:,1) == kk & stimulusConditions(:,2) == jj;
                detectionSuccess(:,kk,jj,ii) = success(idx);
            end
        end
    end
    
    %%
    
    detectionPerf = 100*squeeze(mean(detectionSuccess));
    
    %%
    
    save(resultsFile,'detectionPerf','detectionSuccess');
    
    %%
    
    load(resultsFile);
    
    %%
    
    figure;
    errorbar(repmat(xvals,1,2),squeeze(mean(detectionPerf)),squeeze(std(detectionPerf)));
    line(xlims,[50 50],'Color','k','LineStyle','--');
    line(xlims,[75 75],'Color','k','LineStyle',':');
    
    legend({'Control' 'Drug' 'Chance' 'Threshold'},'Location',legendLocation);
    
    %%
    
    xlabel(xlab);
    xlim(xlims);
    
    ylabel('Decoder Performance (%)');
    ylim([25 100])
    
    set(gca,'XScale',xscale,'XTick',xticks,'XTickLabel',xticklabels);
    
    %%
    
    figFile = ['gratings_perf_vs_' independentVariable '_decoder_' decoder saveFileSuffix];
    saveas(gcf,figFile,'fig');
    saveas(gcf,figFile,'png');
    close(gcf);
end