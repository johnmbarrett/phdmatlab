function combinePSDs(recordings,varargin)
    nRecordings = numel(recordings);
    
    for ii = 1:nRecordings
        [filedir,filename] = getAnalysisOutputDir(recordings{ii});
        
        load(sprintf('%s\\%s_psd.mat',filedir,filename),'psds','ts','F','electrodes');
        
        if ii > 1
            assert(nf == numel(F),'Number of frequencies must not change');
            assert(nt == numel(ts),'Number of time bins must not change');
            assert(nChannels == size(psds,3),'Number of channels must not change');
        else
            psdss = zeros([size(psds) nRecordings]);
        end
        
        nf = numel(F);
        nt = numel(ts);
        nChannels = size(psds,3);
        psdss(:,:,:,ii) = psds;
    end
        
    figure;
    set(gcf,'Renderer','zbuffer','Position',[0 0 1600 900]);
    
    [rows,cols] = subplots(nRecordings);
    
    options = getopt('normaliseallchannels=false smoothing=none sdt=1 sdf=1',varargin{:});
    
    if strcmpi(options.smoothing,'gauss')
        freqRes = median(diff(F));
        sdt = options.sdt/freqRes;
        sdf = options.sdf/freqRes;
        
        [Y,X] = ndgrid(-5*sdf:sdf*5,-5*sdt:sdt*5);
        kernel = gauss2d(X,Y,sdt,sdf,0);
        kernel = kernel/sum(sum(kernel));
        
        for jj = 1:nChannels
            for ii = 1:nRecordings
                psdss(:,:,jj,ii) = conv2(psdss(:,:,jj,ii),kernel,'same');
            end
        end
    end
    
    if all(logical(options.normaliseallchannels))
        maxPs = max(max(max(max(psdss))));
    else
        maxPs = squeeze(max(max(max(psdss,[],1),[],2),[],4));
    end
    
    [~,currentDir] = fileparts(pwd);
        
    for jj = 1:nChannels
        clf;
        
        for ii = 1:nRecordings
            subplot(rows,cols,ii);
            psd = squeeze(psdss(:,:,jj,ii));
            psd = conv2(psd,kernel,'same');
            surf(ts,F,psd);
            view(2);
            shading interp;
            xlim(ts([1 nt]));
            ylim(F([1 nf]));
            xlabel('Time/s');
            ylabel('Frequency/Hz');
            title(recordings{ii});
            caxis([0 maxPs(jj)]);
        end
        
        channelLabel = electrodes(jj).EntityLabel(end-1:end);
        figfile = sprintf('%s_channel_%s_psds_smooth_%s_sdt_%f_sdf_%f',currentDir,channelLabel,options.smoothing,options.sdt,options.sdf);
        saveas(gcf,sprintf('%s.fig',figfile));
        saveas(gcf,sprintf('%s.png',figfile));
    end
end