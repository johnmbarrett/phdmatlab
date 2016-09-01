function plotSpontaneousSpectrograms(dir)
    if nargin < 1
        dirs = {};
    elseif iscell(dir)
        dirs = dir;
    else
        dirs = {dir};
    end

    while isempty(dirs)
        dirs = uipickfiles('Prompt','Choose folders containing analysed data');
    end
    
    if ~iscell(dirs)
        return;
    end
    
    currentDir = pwd;
    
    for jj = 1:numel(dirs)
        dir = dirs{jj};
        
        try
            cd(dir);
        catch error
            fprintf('Change to directory %s failed with error %s\n',dir,error.message);
            continue;
        end

        slashes = strfind(dir,'\');
        lastdirname = dir(slashes(end)+1:end);

        load([regexprep(lastdirname,'^[0-9]+_','') '_metaInfo.mat'],'channelInfo');
        load('SNR.mat');

        psdfig = figure;
        specfig = figure;

        for ii = 1:60
            channel = str2num(channelInfo(ii).label);
            plotnum = floor(channel/10)+8*(mod(channel,10)-1);

            figure(psdfig);
            subplot(8,8,plotnum);

            psd = psds{ii};
            pmax = max(max(psd));
            pmin = min(min(psd));

            if ~any(isempty([pmax pmin])) && ~any(isnan([pmax pmin])) && pmax ~= pmin
                cdata = 255*(psd-pmin)/(pmax-pmin);
                image(spectrogramFreqs{ii},spectrogramTimes{ii},cdata');
            end
            
            figure(specfig);
            subplot(8,8,plotnum);
            
            spectrum = spectra{ii};
            freqs = frequencies{ii};
            
            if isempty(spectrum)
                continue;
            end
            
            maxIndex = find(freqs <= 300,1,'last');
            plot(freqs(1:maxIndex),spectrum(1:maxIndex));
            xlim([0 300]);
            ylim([0 max(spectrum(2:maxIndex))]); % ignore massive DC offset if it exists
        end
    end
    
    cd(currentDir);
end