function [lambdas,nCells] = getFiringStatistics(spontRecording,varargin)
    [spontDir,spontFile] = getAnalysisOutputDir(spontRecording);
    
    statsFile = sprintf('%s\\%s_isi_stats.mat',spontDir,spontFile);
    
    options = getopt('overwrite=''no'' maxt=Inf',varargin{:});
    if ~strcmp(options.overwrite,'yes') && exist(statsFile,'file')
        load(statsFile);
        return;
    end
    
%     N = (timeSpan-0.5)*sampleRate;
%     t = randi(floor(N),[n 1])/sampleRate;
    lambdas = cell(60,1);
%     lambdas = cell(60,1);
    nCells = 0;
    
    figure;
    function fn(spikeTimes,channelIndex,channelLabel,cluster,~,varargin)
        nCells = nCells+1;
        
        if numel(spikeTimes) < 2
            if numel(spikeTimes) == 0
                if ~isfinite(options.maxt)
                    % can't do anything
                    lambdas{channelIndex}(cluster) = 0;
                    return;
                end
            
                % find lambda such that the likelihood of lambda given no
                % events in [0,maxT] is 95%
                lambdas{channelIndex}(cluster) = -log(0.95)/maxT;
                return;
            end
            
            if ~isfinite(options.maxt)
                % the most we know is there are no spikes within
                % [0,spikeTimes(1))
                lambdas{channelIndex}(cluster) = -log(0.95)/spikeTimes(1);
                return;
            end
            
            T0 = spikeTimes(1);
            T1 = options.maxt-spikeTimes(1);
            T = max(T0,T1);
            
            % There is no real-valued solution to the problem of finding
            % lambda such that the probability of one spike in any given
            % interval is 95%.  However, the likelihood of lambda given no
            % spikes in T0 or T1 is the probability of no spikes in T0 or 
            % T1 given lambda.  Because interevent intervals are indepen-
            % dent, this probability is the product of the probability of 
            % no spikes in T0 and the probability of no spike in T1.  The 
            % larger the interval, the lower the probability.  Hence, if we
            % ensure the probability of no spikes in the larger interval is
            % sqrt(0.95), then the likelihood >= sqrt(0.95)*sqrt(0.95) = 0.95.
            % That's a good enough lower bound.
            lambdas{channelIndex}(cluster) = -log(sqrt(0.95))/T;
            return;
        end
        
        isis = diff(spikeTimes);
        
        clf;
        [f,x] = hist(isis,sqrt(numel(isis)));
        bar(x,f/sum(f));
        xlabel('ISI/s');
        ylabel('Frequency');
        title(sprintf('%s - Channel %s cluster %d',spontFile,channelLabel,cluster));
        filename = sprintf('%s\\isi_hist_%s_channel_%s_cluster_%d',spontFile,spontFile,channelLabel,cluster);
        saveas(gcf,filename,'fig');
        saveas(gcf,filename,'png');
        
        clf;
        scatter(isis(1:end-1),isis(2:end));
        xlabel('ISI{_i}/s');
        ylabel('ISI{_(i+1)}/s');
        title(sprintf('%s - Channel %s cluster %d',spontFile,channelLabel,cluster));
        filename = sprintf('%s\\isi_correlation_%s_channel_%s_cluster_%d',spontFile,spontFile,channelLabel,cluster);
        saveas(gcf,filename,'fig');
        saveas(gcf,filename,'png');
        
        clf;
        hold on;
        
        cdfplot(isis);
        
        mu = expfit(isis);
        lambda = 1/mu;
        
        % this made sense when I was writing it for analysing ULED
        % responses, but in the general case it should be up to the calling
        % function how to treat extremely low firing rates
%         if lambda < 1
%             lambdas{channelIndex}(cluster) = 0;
%         else
            lambdas{channelIndex}(cluster) = lambda;
%         end
            
        fplot(@(x) expcdf(x,mu),xlim,'Color','r');
        
        phat = gamfit(isis);
        fplot(@(x) gamcdf(x,phat(1),phat(2)),xlim,'Color','g');
        
        parmhat = lognfit(isis);
        fplot(@(x) logncdf(x,parmhat(1),parmhat(2)),xlim,'Color','m');
        
        legend({ ...
            'Empirical CDF' ...
            sprintf('Fitted exponential CDF (%s = %2.3f)','{\lambda}',lambda) ...
            sprintf('Fitted Gamma CDF (a = %2.3f, b = %2.3f)',phat(1),phat(2)) ...
            sprintf('Fitted log-normal CDF (%s = %2.3f, %s = %2.3f)','{\mu}',parmhat(1),'{\sigma}',parmhat(2)) ...
            },'Location','SouthEast');
        title(sprintf('%s - Channel %s cluster %d',spontFile,channelLabel,cluster));
        
        filename = sprintf('%s\\isi_distribution_%s_channel_%s_cluster_%d',spontFile,spontFile,channelLabel,cluster);
        saveas(gcf,filename,'fig');
        saveas(gcf,filename,'png');
%         [f,x] = hist(isis,sqrt(200));
%         bar(x,f);
%         maxF = max(f);
%         fplot(@(x) maxF*exppdf(x,muhat)/exppdf(0,muhat),xlim,'Color','r');
%         if isempty(samples{channelIndex})
%             samples{channelIndex} = zeros(n/10,maxCluster);
%         end
%         
        
%         sample = zeros(n/10,1);
%         
%         for ii = 1:n
%             spikes = spikeTimes(spikeTimes > t(ii) & spikeTimes <= t(ii)+0.5);
%             sample(ceil(ii/10)) = sample(ceil(ii/10)) + numel(spikes)/10;
%         end
%         
%         samples{channelIndex}(:,cluster) = sample;
%         lambdas{channelIndex}(cluster) = poissfit(sample);
    end

    forEachChannel(spontFile,[],true,@fn);
    
    save(statsFile,'lambdas','nCells');
end