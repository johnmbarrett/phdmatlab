function detectGratingResponsiveCells(recording,varargin)
    [fileDir,filename] = getAnalysisOutputDir(recording);
    
    load(sprintf('%s\\%s_psrh_abs.mat',fileDir,filename));
    load(sprintf('%s\\%s_spikes_concat_forceclustered_no_ignorenoise_no.mat',fileDir,filename))
    %%
    
    cells = [str2num(char(cells(:,1:2))) cells(:,3)]; %#ok<NODEF,ST2NM>
    units = cells;

    [~,si] = sort(valuess{2}); %#ok<USENS>
    rasters = squeeze(rasters(cells(:,2) ~= 0,valuess{1} == 64,si,1)); %#ok<NODEF>
    cells = cells(cells(:,2) ~= 0,:);
    nCells = size(cells,1);

%     firingRate = zeros(50,8,2,nCells);
    firingRate = zeros(25,8,nCells);

    for ii = 1:nCells
        for jj = 1:8
            raster = rasters{ii,jj};

%             for kk = 1:3:148
            for kk = 1:6:150
                line = raster{kk};

                nSpikes = sum(line > 0 & line <= 0.5);
                firingRate(ceil(kk/6),jj,ii) = nSpikes; %*2;
%                 firingRate(ceil(kk/3),jj,1,ii) = nSpikes*2;

%                 nSpikes = sum(line > 0.5 & line <= 1);
%                 firingRate(ceil(kk/3),jj,2,ii) = nSpikes;
            end
        end
    end
    
    % exclude any cells that have a less than 50% chance of responding to
    % their most preferred stimulus
    lowFiringCells = all(sum(firingRate) < 13);
    firingRate(:,:,lowFiringCells) = 0;

%     % use odd-half to decide which phases to analyse
%     meanFiringRate = squeeze(mean(firingRate(1:2:end-1,:,:,:)));
%     stdFiringRate = squeeze(std(firingRate(1:2:end-1,:,:,:)));
%     
%     figure;
%     set(gcf,'Position',[0 0 1000 500]);
%     titles = {'Grating' 'Mask'};
    
%     for ii = 1:nCells
%         for jj = 1:2
%             subplot(1,2,jj);
%             errorbar(0:0.125:0.875,meanFiringRate(:,jj,ii),stdFiringRate(:,jj,ii)/sqrt(25));
%             title(titles{jj});
%             xlabel('Phase');
%             xlim([0 1]);
%             
%             if jj == 1
%                 ylabel('Firing Rate/Hz');
%             end
%         end
%         
%         suptitle(sprintf('Channel %d cluster %d',cells(ii,1),cells(ii,2)));
%         
%         figFile = sprintf('%s\\frvsphase_%s_channel_%d_cluster_%d',fileDir,filename,cells(ii,1),cells(ii,2));
%         saveas(gcf,figFile,'fig');
%         saveas(gcf,figFile,'png');
%     end
    
%     close(gcf);

    %%

    % ps = zeros(nCells,2);

    % for ii = 1:nCells
    % [~,mn] = min(mfr(:,ii));
    % [~,mx] = max(mfr(:,ii));
    % ps(ii,1) = ranksum(fr(:,mn,ii),fr(:,mod(mn+3,8)+1,ii));
    % ps(ii,2) = ranksum(fr(:,mx,ii),fr(:,mod(mx+3,8)+1,ii));
    % end

    testedCells = find(~lowFiringCells);
    nTests = numel(testedCells);
    ps = zeros(nTests,1);
    sumFiringRate = squeeze(sum(firingRate));
    alpha = pi*(0:7)'/4;

%     antiPhase = mod((1:8)+3,8)+1;

    for ii = 1:nTests
%         dfr = abs(meanFiringRate(:,:,ii) - meanFiringRate(antiPhase,:,ii));
%         [~,mi] = max(dfr(:));
%         [ph,pl] = ind2sub([8 2],mi);
% 
%         % use even half to choose p-values
%         ps(ii) = ranksum(firingRate(2:2:end,ph,pl,ii),firingRate(2:2:end,mod(ph+3,8)+1,pl,ii));
        ps(ii) = circ_rtest(alpha,sumFiringRate(:,testedCells(ii)),pi/4);
    end

    %%
%     g = zeros(nCells,1);
%     g([3 4 9 13 14 16 21 26 34 35 37 41 42 45 48 50 51 52 59 64]) = 1;
    %%
    % disp([(1:nCells)' cells ps < 0.025 ps < 0.005 ps < 0.0005 fdrcorrect(ps) fdrcorrect(min(ps,[],2)) ps < 0.05/(2*nCells) g]);
%     disp([(1:nCells)' cells ps < 0.05 ps < 0.01 ps < 0.001 fdrcorrect(ps) ps < 0.05/(1*nCells) g]);
    responsiveCellIndices = testedCells(fdrcorrect(ps));
    responsiveCells = cells(responsiveCellIndices,:);
    responsiveUnitIndices = zeros(size(responsiveCellIndices));
    
    for ii = 1:size(responsiveCellIndices,1)
        responsiveUnitIndices(ii) = find(ismember(units,responsiveCells(ii,:),'rows'));
    end
    
    saveFile = sprintf('%s\\%s_grating_responsive_cells',fileDir,filename);
    oldFile = [saveFile '.mat'];
    
    if exist(oldFile,'file')
        oldFileProperties = dir(oldFile);
        movefile(oldFile,[saveFile '_' strrep(strrep(oldFileProperties.date,' ','_'),':','-') '.mat']);
    end
    
%     save(saveFile,'cells','units','firingRate','meanFiringRate','stdFiringRate','responsiveCellIndices','responsiveUnitIndices','responsiveCells');
    save(saveFile,'cells','units','firingRate','responsiveCellIndices','responsiveUnitIndices','responsiveCells','testedCells','ps');
end