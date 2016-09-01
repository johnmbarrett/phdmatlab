function sloanLettersAnalysis(ctrlRecording,drugRecording,sequenceFile,squaress,ctrlCellsSuffix,drugCellsSuffix,saveFileSuffix,rastersOnly)
    % channel names increase going right and down
    % image position indices are row-major starting from the bottom, but
    % because of the way the APS projector works the image is flipped
    % upside-down, so the first image is in the top-left of the array, then it
    % moves leftwards first followed by downwards
    
    if nargin < 8
        rastersOnly = false;
    end
    
    if nargin < 7
        saveFileSuffix = '';
    else
        saveFileSuffix = ['_' saveFileSuffix];
    end
    
    if nargin < 6
        drugCellsSuffix = 'drug';
    end
    
    if nargin < 5
        ctrlCellsSuffix = 'ctrl';
    end

    if nargin < 4 || isempty(squaress)
        if exist('./letter squares.mat','file')
            load('./letter squares.mat','squaress');
        else
            squaress = { ...
                [1 1 64 64], ...
                [1 1 64 64], ...
                [1 1 32 32; 33 1 64 21; 1 33 32 64; 33 33 64 64], ...
                [ ...
                    1 1 21 21; 22 1 43 21; 44 1 64 21; ...
                    1 22 21 43; 22 22 43 43; 44 22 64 43; ...
                    1 44 21 64; 22 44 43 64; 44 44 64 64 ...
                ]};
        end
    end

    % colours = distinguishable_colors(9);
    % 
    % figure;
    % 
    % for ii = 3 
    %     clf;
    %     hold on;
    %     I = imread(sprintf('R:\\data\\APSexperiments\\John\\JBOG0054\\Images\\%d.png',6));
    %     I = I(16:16:end,16:16:end,:);
    %     image(I);
    %     
    %     corners = cornerss{ii};
    %     
    %     for jj = 1:size(corners,1)
    %         corner = corners(jj,:);
    %         fill(corner([1 1 3 3]),corner([2 4 4 2]),colours(jj,:),'EdgeColor','none','FaceAlpha',0.25);
    %     end
    %     
    %     xlim([0 64])
    %     ylim([0 64])
    %     
    % %     input('...');
    % end
    % 
    % close(gcf);

    %%

    load('channelNames.mat');
    load('EventNo.mat');
    
    if nargin < 3 || ~exist(sequenceFile,'file')
        load('./letters sequence.mat');
    elseif exist(sequenceFile,'file')
        load(sequenceFile);
    else
        error('Could not find sequence file');
    end
        
    load('spiketimestamps.mat');

    %%

    if nargin < 1
        ctrlRecording = 4;
    end
    
    if nargin < 2
        drugRecording = 7;
    end
    
    recordings = [ctrlRecording drugRecording];

    ctrl = load(['ChR2_cells_letters_' ctrlCellsSuffix '.mat']);
    drug = load(['ChR2_cells_letters_' drugCellsSuffix '.mat']);

    cellIndices = arrayfun(@(S) S.bestUnits,[ctrl drug],'UniformOutput',false);
    nCells = cellfun(@numel,cellIndices);
    
    %%
    
    if rastersOnly
        disp(sum(nCells));
        
        rasterStimuli = ismember(conditionOrder(:,1),[1 2]);
        rasterConditions = conditionOrder(rasterStimuli,1:2); %#ok<*NODEF>
        
        for ii = 1:2
            rasterTimes = EventNo{recordings(ii)}; %#ok<USENS>
            rasterTimes = rasterTimes(4*(find(rasterStimuli)-1)+1);

            for jj = 1:nCells(ii)
                tic;
                figs = rasterPlot(spiketimestamps{cellIndices{ii}(jj)},rasterTimes,letters(rasterConditions(:,2))',[-0.25 0 0.25 0.75],[],true,0.05,rasterConditions(:,1),{'Letter' 'Size'},true); %#ok<USENS>

                for kk = 1:numel(figs)
                    figFile = sprintf('Stim %d\\letters_raster_channel_%s_size_%d',recordings(ii),channelNames{1,cellIndices{ii}(jj)},kk); %#ok<USENS>
                    saveas(figs(kk),figFile,'fig');
                    saveas(figs(kk),figFile,'png');
                end
                toc;
            end
        end
        
        return;
    end

    %%

    channelLabels = cellfun(@(idx) channelNames(1,idx)',cellIndices,'UniformOutput',false);
    channelCoords = arrayfun(@(n) zeros(n,6),nCells,'UniformOutput',false); % x, y, [image index per scale]
    indicesPerPosition = cell(nPositions,2);
    cumPositionsPerScale = [0 cumsum(positionsPerScale)];

    for ii = 1:2
        labels = vertcat(channelLabels{ii}{:});
        x = str2num(labels(:,6:7)); %#ok<ST2NM>
        y = str2num(labels(:,3:4)); %#ok<ST2NM>

        channelCoords{ii}(:,1:2) = [x y];

        for jj = 1:4
            for kk = 1:size(squaress{jj},1)
                inSquare = ...
                    x >= squaress{jj}(kk,1) & y >= squaress{jj}(kk,2) & ...
                    x <= squaress{jj}(kk,3) & y <= squaress{jj}(kk,4);

                assert(all(channelCoords{ii}(inSquare,jj+2) == 0));
                channelCoords{ii}(inSquare,jj+2) = kk;
                indicesPerPosition{cumPositionsPerScale(jj)+kk,ii} = find(inSquare);
            end
        end
    end

    cellsPerPosition = cellfun(@numel,indicesPerPosition);

    %%

%     trialsPerPosition = zeros(nPositions,2);
    stimuli = repmat({zeros(nLetters*nReps,1)},nPositions,2);
    responses = arrayfun(@(n) zeros(nLetters*nReps,n),cellsPerPosition,'UniformOutput',false);

    %%

    scaleOrder = conditionOrder(:,1); %#ok<NODEF>
    isLetter = scaleOrder > 0 & scaleOrder <= nScales;
    scaleOrder = scaleOrder(isLetter);
    nLetterTrials = sum(isLetter);
    letterOrder = conditionOrder(isLetter,2:end);

    %%

    for ii = 1:2
        stimulusTimes = EventNo{recordings(ii)};
        stimulusTimes = stimulusTimes(4*(find(isLetter)-1)+1);
        assert(numel(stimulusTimes) == nLetterTrials);

        trialsPerScale = zeros(1,4);

        for jj = 1:nLetterTrials
            tic;
            scale = scaleOrder(jj);
            trialsPerScale(scale) = trialsPerScale(scale) + 1;
            assert(all(trialsPerScale(scale) <= 250));

            nSpikes = cellfun(@(t) sum(t > stimulusTimes(jj) & t <= stimulusTimes(jj)+0.25),spiketimestamps(cellIndices{ii}));

            for kk = 1:positionsPerScale(scale)
                positionIndex = cumPositionsPerScale(scale)+kk;
                stimuli{positionIndex,ii}(trialsPerScale(scale)) = letterOrder(jj,kk);
                responses{positionIndex,ii}(trialsPerScale(scale),:) = nSpikes(indicesPerPosition{positionIndex,ii});
            end
            toc;
        end
    end

    %%

    for ii = 1:cumPositionsPerScale(end)
        scale = find(cumPositionsPerScale < ii,1,'last');
        position = ii-cumPositionsPerScale(scale);

        for jj = 1:2
            assert(isequal(stimuli{ii,jj},letterOrder(scaleOrder == scale,position)));
        end
    end

    %%

    perfPerPosition = zeros(nPositions,2);

    for ii = 1:2
        for jj = 1:nPositions
            tic;
            s = stimuli{jj,ii};
            R = responses{jj,ii};

            if size(R,2) == 0
                continue;
            end

            perfPerPosition(jj,ii) = simpleBayesianDecoder(s,R,[],[],25,'discrete');
            toc;
        end
    end

    %%

    datas = {cellsPerPosition 100*perfPerPosition};
    xlabels = {'# Cells' 'Performance (%)'};
    ylabels = {'Control' 'Drug'};

    for hh = 1:2
        figure;
        set(gcf,'Position',[100 100 800 400]);

        for jj = 1:4
            data = datas{hh}(cumPositionsPerScale(jj)+1:cumPositionsPerScale(jj+1),:);

            for ii = 1:2
                subplot('Position',[0.05+0.24*(jj-1) 0.025+0.55*(2-ii) 0.2 0.35+0.15*(ii-1)]);    
                rows = rowsPerScale(jj);
                surf([flipud(reshape(data(:,ii),rows,rows)') zeros(rows,1); zeros(1,rows+1)]);
                caxis([0 max(max(data))]);
                set(gca,'XTick',[],'YTick',[]);
                view(2);

                if ii == 1
                    title(sprintf('%dx%d letters',rows,rows));
                end

                if ii == 2
                    xlabel(xlabels{hh});
                end

                if jj == 1
                    ylabel(ylabels{ii});
                end

                xlim([0 rows]+1);
                ylim([0 rows]+1);
            end

            colorbar('SouthOutside');
        end
    end

    %%

    saveas(gcf,['letter_perf_vs_pos' saveFileSuffix],'fig');
    saveas(gcf,['letter_perf_vs_pos' saveFileSuffix],'png');
    close(gcf);

    saveas(gcf,['letter_ncells_vs_pos' saveFileSuffix],'fig');
    saveas(gcf,['letter_ncells_vs_pos' saveFileSuffix],'png');
    close(gcf);

    %%

    figure;
    scatter(cellsPerPosition(:),perfPerPosition(:))
    xlabel('# Cells');
    ylabel('Performance (%)');
    set(gca,'XScale','log','YScale','log');

    %%

    saveas(gcf,['letter_ncells_vs_perf' saveFileSuffix],'fig');
    saveas(gcf,['letter_ncells_vs_perf' saveFileSuffix],'png');
    close(gcf);
    
    %%
    
    perfPerScale = zeros(nScales,2);
    
    for ii = 1:nScales
        perfPerScale(ii,:) = mean(perfPerPosition(cumPositionsPerScale(ii)+1:cumPositionsPerScale(ii+1),:),1);
    end
    
    %%
    
    figure;
    bar(100*perfPerScale)
    ylim([0 100])
    line(xlim,[10 10],'Color','k','LineStyle','--');
    legend({'Control' 'Drug' 'Chance'},'Location','NorthEast')
    letterVisualAngles = sizes*4/30;
    set(gca,'XTickLabel',letterVisualAngles)
    xlabel('Letter Size (degress of visual angle)')
    ylabel('Decoder Performance (%)')

    %%

    saveas(gcf,['letter_perf_vs_scale' saveFileSuffix],'fig');
    saveas(gcf,['letter_perf_vs_scale' saveFileSuffix],'png');
    close(gcf);

    %%

    [rho,pval] = corr(cellsPerPosition(:),perfPerPosition(:)); %#ok<NASGU,ASGLU>

    %%

    save(['letter_responses' saveFileSuffix '.mat'],'cellIndices','cellsPerPosition','channelCoords','channelLabels','indicesPerPosition','letterOrder','perfPerPosition','pval','responses','rho','scaleOrder','squaress','stimuli');
end