function analyseGratingInfo(stimFile,excludeTrials,recomputeInfoPerf,recomputePID)
    if nargin < 4
        recomputePID = true;
    end
    
    if nargin < 3
        recomputeInfoPerf = true;
    end
   
    if nargin < 2
        excludeTrials = (1:6:150)';
    end
    
    if nargin < 1
        stimFile = 'Stim 2';
    end
    
    %%

    load(sprintf('./%s_stp_data.mat',stimFile));

    %%
    
    nCells = size(responsiveCells,1);
    
    if recomputePID || recomputeInfoPerf
        counts = cellfun(@numel,responses); %#ok<NODEF>
        latencies = arrayfun(@(t,c) tertiaryfun(c,@() t{1}(1),@() Inf),responses,counts); %[]),responses,counts,'UniformOutput',false);
        linf = arrayfun(@(t,c) tertiaryfun(c,@() t{1}(1),@() []),responses,counts,'UniformOutput',false);
    end

    %%

    % polarities = {'On' 'Off'};
    % 
    % for ii = 1:2
    %     for jj = 1:4
    %         figure;
    %         for kk = 1:nCells
    %             subplot(2,4,kk);
    %             ch = arrayfun(@(s) median(counts(stimuli(:,jj) == s,ii,jj,kk)),unique(stimuli(:,jj)));
    %             lh = arrayfun(@(s) median(latencies(stimuli(:,jj) == s,ii,jj,kk)),unique(stimuli(:,jj)));
    %             bar([ch lh]);
    %         end
    %         
    %         suptitle(sprintf('%s Responses, Confounder = %d',polarities{ii},jj));
    %     end
    % end

    %%

    subunits = [2.^(0:nextpow2(nCells)-1) nCells];
    saveFile = sprintf('%s_info_results.mat',stimFile);
    
    if exist(saveFile,'file')
        load(saveFile);
    else
        recomputeInfoPerf = true;
        recomputePID = true;
        infos = zeros(nCells,4,2,4); % cell, confounding variable, on/off, coding strategy
        perfs = zeros(10,numel(subunits),4,2,5);
    end
    
    %%
    
    n = size(responses,4);
    disp(n);
    k = size(responses,1)/8;
    excludeStimuli = repmat(excludeTrials,8,1) + kron(k*(0:7)',ones(numel(excludeTrials),1));
    includeStimuli = setdiff((1:size(responses,1))',excludeStimuli);
    N = n*(n-1)/2;
    disp(N);
    
    %%
    
    if true
        ws1 = warning('off','GetOpt:UnknownOption'); 

        dI = zeros(N,2,4);

        for hh = 1:4
            nn = 0;
            for ii = 1:n-1
                for jj = ii+1:n
                    tic;
                    nn = nn + 1;
                    dI(nn,1,hh) = biasCorrect(stimuli(includeStimuli,hh),squeeze(counts(includeStimuli,1,hh,[ii jj])),'method','function','mifunction',@deltaI);
                    dI(nn,2,hh) = biasCorrect(stimuli(includeStimuli,hh),squeeze(counts(includeStimuli,1,hh,[ii jj])));
                    toc;
                end
            end
        end

        save(saveFile,'-append','dI');

        warning(ws1);

%         retursn;
    end
    
    %%
    
%     load(saveFile);
    
    if ~exist('pidt','var') || recomputePID
        info1 = zeros(n,4);
        kld1 = zeros(n,8,4);

        for ii = 1:n
            for jj = 1:4
                tic;
                r = squeeze(responses(includeStimuli,1,jj,ii));
                kld1(ii,:,jj) = binlessInfo(stimuli(includeStimuli,jj),r,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1,'kld',true); %#ok<NODEF>
                info1(ii,jj) = binlessInfo(stimuli(includeStimuli,jj),r,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1);
                toc;
            end
        end

        info1 = max(info1,0);

        info2 = zeros(N,4);
        pairs = zeros(N,2);

        nn = 0;
        for ii = 1:n-1
            for jj = ii+1:n
                nn = nn + 1;

                pairs(nn,:) = [ii jj];

                for kk = 1:4
                    tic;
                    r = squeeze(responses(includeStimuli,1,kk,[ii jj]));
                    info2(nn,kk) = binlessInfo(stimuli(includeStimuli,kk),r,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1);
                    toc;
                end
            end
        end

        info2 = max(info2,0);
    else
        pairs = zeros(N,2);

        nn = 0;
        for ii = 1:n-1
            for jj = ii+1:n
                nn = nn + 1;

                pairs(nn,:) = [ii jj];
            end
        end
    end
    
    red = squeeze(mean(min(cat(4,kld1(pairs(:,1),:,:),kld1(pairs(:,2),:,:)),[],4),2));
    uq1 = info1(pairs(:,1),:,1)-red;
    uq2 = info1(pairs(:,2),:,1)-red;
    syn = info2-uq1-uq2-red;
    pidt = permute(cat(3,red,uq1,uq2,syn,info2),[1 3 2]); %#ok<NASGU>
    
    save(saveFile,'-append','pidt','info1','info2','kld1');
    
    return;
    
    %%
    
    n = size(responses,1);
    
    %%
    
    if recomputeInfoPerf
        [trainDecoder,initDecoder] = getJacobsDecoderTrainingFun(0.001,0.025);
        jacobsDecoder = struct('initDecoder',initDecoder,'trainDecoder',trainDecoder,'likelihoodFun',@getJacobsLikelihood);
        wfsDecoder = struct('trainDecoder',@trainWFSDecoder,'makeSubDecoder',@makeWFSSubDecoder,'likelihoodFun',@getWFSLikelihood);

        ws = warning('query','GetOpt:UnknownOption');
        warning('off','GetOpt:UnknownOption');

        for ii = 1:2
            for jj = 1:4
    %             tic;
                [~,~,x] = unique(stimuli(includeStimuli,jj)); %#ok<NODEF>
    %             
    %             for kk = 1:nCells
    %                 [~,~,y] = unique(counts(includeStimuli,ii,jj,kk));
    %                 infos(kk,jj,ii,1) = discreteMutualInformation(x,y)/log(2); %biasCorrect(stimuli(:,jj),counts(:,ii,jj,kk));
    %                 infos(kk,jj,ii,2) = discreteMutualInformation(x,(counts(includeStimuli,ii,jj,kk) > 0)+1)/log(2); %biasCorrect(stimuli(:,jj),counts(:,ii,jj,kk));
    %                 infos(kk,jj,ii,3) = binlessInfo(stimuli(includeStimuli,jj),linf(includeStimuli,ii,jj,kk),'variabletype','mixed','triallength',max(stimulusLengths(jj,:,ii)),'label',sprintf('Ch%d_%d',responsiveCells(kk,1),responsiveCells(kk,2)),'samplingfrequency',25000,'stratstrat',1,'singstrat',0,'maxembed',1);
    %                 infos(kk,jj,ii,4) = binlessInfo(stimuli(includeStimuli,jj),responses(includeStimuli,ii,jj,kk),'variabletype','mixed','triallength',max(stimulusLengths(jj,:,ii)),'label',sprintf('Ch%d_%d',responsiveCells(kk,1),responsiveCells(kk,2)),'samplingfrequency',25000,'stratstrat',1,'singstrat',0,'maxembed',max(counts(:)));
    %             end
    %             toc;
    % 
    %             tic;
    %             [~,~,perfs(:,:,jj,ii,1)] = simpleBayesianDecoder(x,permute(counts(includeStimuli,ii,jj,:),[1 4 2 3]),[],subunits,10,'discrete');
    %             toc;
    %             tic;
    %             [~,~,perfs(:,:,jj,ii,2)] = simpleBayesianDecoder(x,(permute(counts(includeStimuli,ii,jj,:),[1 4 2 3]) > 0)+1,[],subunits,10,'discrete');
    %             toc;
    %             tic;
    %             [~,~,perfs(:,:,jj,ii,3)] = simpleBayesianDecoder(x,permute(latencies(includeStimuli,ii,jj,:),[1 4 2 3]),[],subunits,10,'mixedBernoulliKDE',0,max(stimulusLengths(jj,:,ii)));
    %             toc;
    %             tic;
    %             [~,~,perfs(:,:,jj,ii,4)] = simpleBayesianDecoder(x,permute(responses(includeStimuli,ii,jj,:),[1 4 2 3]),[],subunits,10,jacobsDecoder);
    %             toc;
                tic;
                [~,~,perfs(:,:,jj,ii,5)] = simpleBayesianDecoder(x,permute(latencies(includeStimuli,ii,jj,:),[1 4 2 3]),[],subunits,10,wfsDecoder);
                toc;
            end
        end

        warning(ws);
    end

    %%

    polarities = {'On' 'Off'};
    codes = {'Count' 'Response' 'Latency' 'Timing' 'WFS'};

    minfo = squeeze(median(infos,1));
    figure
    set(gcf,'Position',[100 100 800 400]);
    for ii = 1:2
        subplot(1,2,ii);
        bar(squeeze(minfo(:,ii,:))');
        set(gca,'XTickLabel',codes(1:4));
        title(polarities{ii});
        ylabel('Mutual Information (bits)');
    end

    %%

    figFile = sprintf('%s_info_v_code',stimFile);
    saveas(gcf,figFile,'fig');
    export_fig(figFile,'-png','-transparent','-painters','-m1');
    close(gcf);

    %%

    nperfs = 100*permute(mean(perfs,1),[2:5 1]);

    yy = [5*floor(min(nperfs(:))/5) 5*ceil(max(nperfs(:))/5)];

    figure;
    set(gcf,'Position',[100 100 1200 600]);
    for jj = 1:2
        for ii = 1:4
            subplot(2,4,4*(jj-1)+ii);
            plot(subunits,nperfs(:,:,jj,ii));
            line(xlim,[12.5 12.5],'Color','k','LineStyle','--');
            
            set(gca,'XScale','log');

            title([codes{ii} ' - ' polarities{jj}]);

            xlabel('# Neurons');

            if ii == 1
                ylabel('Decoder performance (%)');
            end

            ylim(yy);
        end
    end

    %%

    figFile = sprintf('%s_perf_v_neurons',stimFile);
    saveas(gcf,figFile,'fig');
    export_fig(figFile,'-png','-transparent','-painters','-m1');
    close(gcf);

    %%

    mperfs = permute(nperfs(end,:,:,:),[2:4 1]);

    figure
    set(gcf,'Position',[100 100 800 400]);
    for ii = 1:2
        subplot(1,2,ii);
        bar(squeeze(mperfs(:,ii,:))');
        set(gca,'XTickLabel',codes);
        title(polarities{ii});
        ylabel('Decoder Performance (%)');
    end

    %%

    figFile = sprintf('%s_perf_v_code',stimFile);
    saveas(gcf,figFile,'fig');
    export_fig(figFile,'-png','-transparent','-painters','-m1');
    close(gcf);

    %% 

    if recomputePID
        pids = pairwisePID(stimuli(:,1),permute(counts,[4 1:3]));
    end

    %%

    mpid = permute(median(pids,1),[2 4 3 1]);

    figure
    set(gcf,'Position',[100 100 1200 500]);

    for ii = 1:2
        subplot(1,2,ii);
        bar(squeeze(mpid(:,:,ii)));
        set(gca,'XTickLabel',{'Redundant' 'Unique 1' 'Unique 2' 'Synergistic' 'Total Info'});
        title(polarities{ii});
        ylabel('Mutual Information (bits)');
    end

    %%

    figFile = sprintf('%s_pid',stimFile);
    saveas(gcf,figFile,'fig');
    export_fig(figFile,'-png','-transparent','-painters','-m1');
    close(gcf);

    %%

    save(saveFile,'-append','counts','infos','latencies','linf','perfs','pids','subunits');
end