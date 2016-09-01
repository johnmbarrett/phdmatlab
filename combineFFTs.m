function combineFFTs(recordings,varargin)
    options = getopt('window=''none'' points=0 suffix=NaN',varargin{:});

    for ii = 1:numel(recordings)
        recording = recordings{ii};
        
        [filedir,filename] = getAnalysisOutputDir(recording);
        
        fftfile = sprintf('%s\\%s_fft',filedir,filename);
        load(fftfile,'specs','freqs','interesting','electrodes');
        
        if ii == 1
            allSpecs = zeros(size(specs,1),numel(recordings),size(specs,3));
            filteredSpecs = zeros(size(specs,1),numel(recordings),size(specs,3));
            allFreqs = freqs(interesting);
        else
            assert(isequal(allFreqs,freqs(interesting)));
        end
        
        allSpecs(:,ii,:) = specs;
        
        windowFn = sprintf('%swin',options.window);
        if exist(windowFn,'file')
            eval(sprintf('win = %s(%d);',windowFn,options.points));
            win = win/sum(win);
            
            for jj = 1:size(specs,3)
                specs(:,1,jj) = conv(specs(:,1,jj),win,'same'); %#ok<AGROW>
            end
        end
        
        filteredSpecs(:,ii,:) = specs;
    end
    
    [~,currentDir] = fileparts(pwd);
    suffix = currentDir;
    
    if ischar(options.suffix)
        suffix = [suffix '_' options.suffix];
    end
    
    peakFreqs = zeros(size(allSpecs,3),size(allSpecs,2));
    controlPowers = zeros(size(allSpecs,3),size(allSpecs,2));
    peakPowers = zeros(size(allSpecs,3),size(allSpecs,2));
    
    specFigure = figure;
    powerFigure = figure;
    for ii = 1:size(allSpecs,3)
        spec = squeeze(allSpecs(:,:,ii));
        
        figure(specFigure);
        clf;
        hold on;
        colourOrder = distinguishable_colors(size(spec,2));
        set(gca,'ColorOrder',colourOrder);
        plot(allFreqs,filteredSpecs(:,:,ii));
        legend(recordings);
        xlim([0 45]);
        xlabel('Frequency/Hz');
        ylabel('Amplitude/V');
        channelLabel = electrodes(ii).EntityLabel(end-1:end);
        title(sprintf('Channel %s Frequency Spectrum',channelLabel));
        
        % originally only considered data with 'significant' oscillations,
        % but this is double-dipping, so from now on always take all
        % channels that were analysed
%         Z = zscore(spec);
        [maxSpec,imax] = max(spec);
        
%         if any(maxZ >= 6)
            for jj = 1:numel(maxSpec)
                index = imax(jj);
                plot(allFreqs(index),filteredSpecs(index,jj,ii),'Color',colourOrder(jj,:),'LineStyle','none','Marker','o');
                peakFreqs(ii,jj) = allFreqs(index);
                peakPowers(ii,jj) = getPeakPower(spec(:,jj),imax(jj));
            end
        
            figure(powerFigure);
            clf;
            hold on;
        
            plot(peakPowers(ii,:),'Color','b');
            set(gca,'XTick',1:numel(recordings),'XTickLabel',recordings);
            rotateXLabels(gca,45);
            ylabel('Power');
            title(sprintf('Channel %s Oscillation Strength',channelLabel));
            legend('Peak within recording');
%         end
        
        figfile = sprintf('%s_channel_%s_fft_window_%s_points_%d',suffix,channelLabel,options.window,options.points);
        saveas(specFigure,figfile,'fig');
        saveas(specFigure,figfile,'png');
                
%         if maxZ(1) >= 6
            controlPowers(ii,:) = getPeakPower(spec,imax(1));
            plot(controlPowers(ii,:),'Color','r')
            legend({'Peak within recording' 'Control peak'});
%         end
        
        figfile = sprintf('%s_channel_%s_peak_power',suffix,channelLabel,options.window,options.points);
        saveas(powerFigure,figfile,'fig');
        saveas(powerFigure,figfile,'png');
    end
    
    close(specFigure);
    close(powerFigure);
    
    figure;
    boxplot(controlPowers,'labels',recordings,'labelorientation','inline');
    ylabel('Power');
    title(sprintf('Oscillation Strength at Peak Frequency During Control'));
    
    figfile = sprintf('%s_control_power',suffix);
    saveas(gcf,figfile,'fig');
    saveas(gcf,figfile,'png');
    
    boxplot(peakPowers,'labels',recordings,'labelorientation','inline');
    title(sprintf('Oscillation Strength at Peak Frequency in Each Recording'));
    
    figfile = sprintf('%s_peak_power',suffix);
    saveas(gcf,figfile,'fig');
    saveas(gcf,figfile,'png');
    
    xlsxfile = sprintf('%s_fft_data.xlsx',suffix);
    
    if exist(xlsxfile,'file')
        delete(xlsxfile);
    end
    
    worksheets = {'Peak Frequency' 'Peak Power (control peak)' 'Peak Power (local peak)'};
    % powers are in V (or V/s or something) but are on the order of uV, so
    % get them into a range where Excel is unlikely to have underflow probs
    datas = {peakFreqs controlPowers*1e6 peakPowers*1e6};
    
    for ii = 1:numel(worksheets)
        lastColumn = char('A'+numel(recordings)-1);
        xlswrite(xlsxfile,recordings(:)',worksheets{ii},sprintf('A1:%s1',lastColumn));
        xlswrite(xlsxfile,datas{ii},worksheets{ii},sprintf('A2:%s%d',lastColumn,size(datas{ii},1)+1));
    end
end

function power = getPeakPower(spec,index)
    peakStart = index-find(spec(index-1:-1:1,1) < spec(index,1)/2,1);
    peakEnd = find(spec(index+1:end,1) < spec(index,1)/2,1)+index;
    power = sum(spec(peakStart:peakEnd,:));
end