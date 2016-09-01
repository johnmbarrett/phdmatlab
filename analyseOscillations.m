function analyseOscillations(recording,varargin)
    if strcmp(recording(end-3:end),'.raw')
        file = openMCSRawFile(recording);
        sampleRate = file.sampleRate;
        nChannels = numel(file.electrodes);
        nSamples = file.nSamples;
        electrodes = file.electrodes; %#ok<NASGU>
        isRaw = true;
    else
        [file,fileInfo] = openMCDFile(recording);
        sampleRate = round(1/fileInfo.TimeStampResolution);
        [channelIndices,channels] = getMCDChannels(file,fileInfo,'el',varargin{:});
        channelIndices = channelIndices{1};
        channels = channels{1};
        nChannels = numel(channelIndices);
        if nChannels < 1
            disp('Nothing to do');
            return;
        end
        nSamples = channels(1).ItemCount;
        electrodes = channels; %#ok<NASGU>
        isRaw = false;
    end
%     destructor = onCleanup(@() fclose(file.handle));
    
    
    opt = getopt('chunk=60 dofft=no fftsamples=2^21 fftmin=0.1 fftmax=45 fftthresh=6 timestep=5 dospec=no maxlfpfreq=150 window=2^17 freqs=1:30 ytick=0:5:30 overwrite=no',varargin);
    
    [fileDir,filename] = getAnalysisOutputDir(recording);
    
    if strcmpi(opt.dofft,'yes')
        savefile = sprintf('%s\\%s_fft.mat',fileDir,filename);
        N = opt.fftsamples;
        
        tTotal = 0;
        if strcmpi(opt.overwrite,'yes') || ~exist(savefile,'file')
            if isRaw
                nSpecs = floor(nSamples/N);
            else
                nSpecs = 1;
            end

            peakFreqs = cell(nSpecs,1);

            freqs = sampleRate/2*linspace(0,1,N/2+1);
            interesting = freqs >= opt.fftmin & freqs <= opt.fftmax;
            interestingFreqs = freqs(interesting);

            nFreqs = sum(interesting);
            specs = zeros(nFreqs,nSpecs,nChannels);

            iSpecs = 1;
            iChannel = 1;
            while true
                tic;
                if isRaw
                    [data,~,eof] = getMCSRawData(file,1,N,[],'samples',0);
                else
                    channelIndex = channelIndices(iChannel);
                    if N > nSamples
                        startIndex = 1;
                        itemCount = nSamples;
                    else
                        startIndex = max(floor((nSamples-N)/2),1);
                        itemCount = N;
                    end
                    [~,~,data] = ns_GetAnalogData(file,channelIndex,startIndex,itemCount);
                    eof = false;
                end
                tData = toc;

                if eof
                    tTotal = tTotal + tData;
                    break;
                end

                fprintf('Fetched chunk %d/%d of file %s in\t\t\t\t%f seconds\n',iSpecs,nSpecs,recording,tData);

                tic;

                spec = fft(data,N)/N;
                spec = 2*abs(spec(1:N/2+1,:));
                spec = spec(interesting,:);

                specs(:,iSpecs,iChannel) = spec;

                tFFT = toc;
                fprintf('Computed FFT for chunk %d/%d of file %s in\t\t%f seconds\n',iSpecs,nSpecs,recording,tFFT);

                tic;

                Z = (spec-repmat(mean(spec),size(spec,1),1))./repmat(std(spec),size(spec,1),1);
                [maxZ,index] = max(Z);

                peakFreqs{iSpecs} = interestingFreqs(index(maxZ > opt.fftthresh));

                tPeak = toc;
                fprintf('Found peak frequency chunk %d/%d of file %s in\t%f seconds\n',iSpecs,nSpecs,recording,tPeak);

                tTotal = tTotal + tData + tFFT + tPeak;
                
                if isRaw
                    iSpecs = iSpecs + 1;
                elseif iChannel == nChannels
                    break;
                else
                    iChannel = iChannel + 1;
                end
            end
            
            tic;
            save(savefile,'peakFreqs','freqs','interesting','specs','electrodes');
        else
            tic;
            load(savefile);
        end
        
        tTotal = tTotal + toc;
        
        tic;
        
        if isRaw
            T = N/sampleRate;

            x = [];
            g = [];
            for ii = 1:numel(peakFreqs)
                x = [x; peakFreqs{ii}(:)]; %#ok<AGROW>
                g = [g; T*ii*ones(numel(peakFreqs{ii}),1)]; %#ok<AGROW>
            end

            figure;
            boxplot(x,g);
            ts = unique(g);
            xlim([0 numel(ts)+1]);
            maxT = ceil(file.nSamples/(60*file.sampleRate));
            set(gca,'XTick',(0:opt.timestep:maxT)*60*file.sampleRate/N,'XTickLabel',0:opt.timestep:maxT);
            xlabel('Time/min');
            ylabel('Peak frequency/Hz');
            title(filename);

            figfile = sprintf('%s\\%s_peakfreq',fileDir,filename);
            saveas(gcf,figfile,'fig');
            saveas(gcf,figfile,'png');
        else
            figure;
            
            for ii = 1:size(specs,3)
                plot(freqs(interesting),squeeze(specs(:,1,ii)));
                xlim([0 45]);
                xlabel('Frequency/Hz');
                ylabel('Amplitude/V');
                channelLabel = channels(ii).EntityLabel(end-1:end);
                title(sprintf('Channel %s Frequency Spectrum',channelLabel));
            
                figfile = sprintf('%s\\%s_channel_%s_fft',fileDir,filename,channelLabel);
                saveas(gcf,figfile,'fig');
                saveas(gcf,figfile,'png');
            end
        end
        
        tTotal = tTotal + toc;
    end
    
    if strcmpi(opt.dospec,'yes')
        savefile = sprintf('%s\\%s_psd.mat',fileDir,filename);

        tTotal = 0;
        if strcmpi(opt.overwrite,'yes') || ~exist(savefile,'file')
            tic;

            if isRaw
                N = nSamples;
            else
                N = opt.fftsamples;
            end
            
            % this seems to work and I don't particularly care why; I think it's
            % something to do with the half-windows at the beginning and end not
            % being there
            nt = floor(N/(opt.window/2))-1;
            psds = zeros(numel(opt.freqs),nt,nChannels);
            ts = zeros(1,nt);
            chunk = opt.chunk * opt.window;
            % TODO : specify filter options
            [b,a] = ellip(2,0.1,40,(2*opt.maxlfpfreq/sampleRate));

            t = 0;
            nChunks = 1;

            tTotal = toc;

            iChannel = 1;
            while true
                tic;
                if isRaw
                    [data,~,eof] = getMCSRawData(file,1,chunk,[],'samples',0);
                else
                    channelIndex = channelIndices(iChannel);
                    if N > nSamples
                        startIndex = 1;
                        itemCount = nSamples;
                    elseif N < nSamples
                        startIndex = max(floor((nSamples-N)/2),1);
                        itemCount = N;
                    end
                    [~,~,data] = ns_GetAnalogData(file,channelIndex,startIndex,itemCount);
                    eof = false;
                end
                tData = toc;
                fprintf('Fetched chunk #%d of file %s in\t\t%f seconds\n',nChunks,recording,tData);

                tic;
                lfp = filtfilt(b,a,data);

                if isRaw && nChunks > 1
                    lfp = [oldlfp(end-opt.window/2+1:end,:); lfp]; %#ok<AGROW>
                end

                oldlfp = lfp;

                tFilter = toc;
                fprintf('Filtered chunk #%d of file %s in\t%f seconds\n',nChunks,recording,tFilter);

                tTotal = tTotal + tData + tFilter;

                for cc = iChannel:(isRaw*nChannels+~isRaw*iChannel)
                    tic;

                    if isRaw
                        chlabel = file.electrodes(cc).label;
                    else
                        chlabel = channels(cc).EntityLabel(end-1:end);
                    end

                    [~,F,T,P] = spectrogram(lfp(:,1+isRaw*(cc-1)),hamming(opt.window),opt.window/2,opt.freqs,sampleRate);
                    assert(isequal(opt.freqs(:),F(:)));
                    psds(:,t+1:t+size(P,2),cc) = P;

                    tSpec = toc;                                                    
                    fprintf('Computed spectrogram chunk #%d for channel %s (%d/%d) in\t%fseconds\n',nChunks,chlabel,cc,nChannels,tSpec);
                    tTotal = tTotal + tSpec;
                end

                tic;
                if isRaw
                    if t > 0
                        ts(t+1:t+size(T,2)) = T+ts(t);
                    else
                        ts(t+1:t+size(T,2)) = T;
                    end

                    t = t+size(P,2);
                else
                    ts = T;
                end

                nChunks = nChunks + 1;
                tTotal = tTotal + toc;
                
                if isRaw
                    if eof
                        break;
                    end
                elseif iChannel == nChannels
                    break;
                else
                    iChannel = iChannel + 1;
                end
            end

            tic;
            save(savefile,'psds','ts','F','opt','electrodes');
        else
            tic;
            load(savefile);
            nChannels = size(psds,3);
        end

        figure;
        set(gcf,'Renderer','zbuffer');

        [X,Y] = meshgrid(-25:25,5:5);
        Z = gauss2d(X,Y,5,1,0,1);

        tTotal = tTotal + toc;

        for cc = 1:nChannels
            tic;

            psd = psds(:,:,cc);

            spsd = conv2(psd,Z,'same');

            clf;
            hold on;

            surf(ts,F,spsd,'EdgeColor','none');

%             cax = [min(psd(isfinite(spsd))) max(psd(isfinite(spsd)))];
%             caxis(cax);

            view(2);

            set(gca,'YTick',opt.ytick);
            xlabel('Time/s');
            ylabel('Frequency/Hz');

            if isRaw
                chlabel = file.electrodes(cc).label;
            else
                chlabel = channels(cc).EntityLabel(end-1:end);
            end
            
            title(sprintf('Spectrogram for Channel %s (%s)',chlabel,filename));

            xlim([0 ts(end)]);
            ylim(opt.freqs([1 end]));

            figfile = sprintf('%s\\%s_channel_%s_psd',fileDir,filename,chlabel);
            saveas(gcf,figfile,'fig');
            saveas(gcf,figfile,'png');

            tPlot = toc;
            fprintf('Plotted data for channel %d/%d in\t%f seconds\n',cc,nChannels,tPlot);

            tTotal = tTotal + tPlot;
        end
    end
    
    fprintf('Finished analysing oscillations for recording %s in %f seconds\n',recording,tTotal);
end