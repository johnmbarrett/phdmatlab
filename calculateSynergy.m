figPrefix = 'P39';
figSuffix = '_bc';
% figSuffix = '_rand_bc';
respFileSuffix = '';
% respFileSuffix = '_rand';
trials = 1:840;
isMCS = false;
isShuffle = false;

%%

if isMCS
%     stimFile = 2;
    load(sprintf('%s_stp_data.mat',stimFile),'stimuli','responseIndices');

    if isShuffle
        firstNeuron = load(sprintf('%s_stp_data.mat',stimFile),'responses');
    else
        firstNeuron = load(sprintf('%s_stp%s_data.mat',stimFile,respFileSuffix),'responses');
    end
    
    secondNeuron = load(sprintf('%s_stp%s_data.mat',stimFile,respFileSuffix),'responses');
    
    n = numel(responseIndices);
    N = n*(n-1)/2;
    indices = 1:n;
    
    stimuli = stimuli*8+1;
    
    szr = size(firstNeuron.responses);
    szr = num2cell(szr);
    resp1 = cellfun(@(R) R(trials,:),squeeze(mat2cell(cellfun(@numel,permute(firstNeuron.responses(:,1,:,:),[1 3 4 2])),szr{[1 3]},ones(n,1))),'UniformOutput',false);
    resp2 = cellfun(@(R) R(trials,:),squeeze(mat2cell(cellfun(@numel,permute(secondNeuron.responses(:,1,:,:),[1 3 4 2])),szr{[1 3]},ones(n,1))),'UniformOutput',false);
else
    load('channelNames.mat')

%%

    load('EventNo.mat');
    load('spiketimestamps.mat');

    % conditions = [kron([128;64;32;16],ones(8,1)) repmat((1:8)',4,1)];
    % disp(conditions);
    % apsToStimulusReponsePairs('srp.mat',spiketimestamps,EventNo(34:39),'info.txt',0,conditions,channelNames);

    %%

    clear idx idxBOTH
    load('ALLunits25.mat','idxBOTH','idx')

    if exist('idxBOTH','var') && ~exist('idx','var')
        idx = idxBOTH;
    end

    if exist('./IDXalternative.mat','file')
        load('./IDXalternative.mat');
        idx = idxGood;
%         figSuffix = '_alt';
    end
    
    rffile = dir(['RF_selected_' figPrefix '*']);
    
    if ~isempty(rffile);
        load(rffile.name,'RF_selected');
        idx = cellfun(@(s) s.id,RF_selected)+1;
    end
    
    %%
    
    load('srp.mat','stimuli');

    if isShuffle
        firstNeuron = load('srp.mat', 'responses');
    else
        firstNeuron = load(sprintf('srp%s.mat',respFileSuffix), 'responses');
    end
    
    secondNeuron = load(sprintf('srp%s.mat',respFileSuffix), 'responses');

    %%

    n = numel(idx);
    N = n/2*(n-1);

    valid = [channelNames{6,:}] == 0;

    assert(all(ismember(idx-1,find(valid))),'Selected units should only included non-marked units');

    indices = zeros(n,1);

    for ii = 1:n
        indices(ii) = find(ismember(find(valid),idx(ii)-1));
    end
    
    resp1 = cellfun(@(R) permute(R(trials,1,1,:),[1 4 2 3]),firstNeuron.responses(indices),'UniformOutput',false);
    resp2 = cellfun(@(R) permute(R(trials,1,1,:),[1 4 2 3]),secondNeuron.responses(indices),'UniformOutput',false);
end

disp(n);
disp(N);

%%

legends = {'75 mcpd' '37.5 mcpd' '18.4 mcpd' '9.4 mcpd'};

% [sources,sets,trans] = PIDLattice(2);

%% NOISE CORRELATIONS

R = zeros(numel(trials)/8,n,8,4);
rho = zeros(n,n,8,4);
prho = zeros(n,n,8,4);
chi2 = zeros(n,n,8,4);
pchi2 = zeros(n,n,8,4);
cv2 = zeros(n,n,8,4);

for ii = 1:4
    for jj = 1:8
        tic;
        for kk = 1:n
            R(:,kk,jj,ii) = resp2{kk}(stimuli(trials,ii) == jj,ii);
        end
        toc;
        
        tic;
        [rho(:,:,jj,ii),prho(:,:,jj,ii)] = corr(R(:,:,jj,ii));
        toc;
        
        for kk = 1:n-1
            for ll = kk+1:n
                tic;
                [t,chi2(kk,ll,jj,ii),pchi2(kk,ll,jj,ii)] = crosstab(R(:,kk,jj,ii),R(:,ll,jj,ii));
                
                cv2(kk,ll,jj,ii) = (chi2(kk,ll,jj,ii)/(numel(trials)/8))/max(1,min(size(t)-1));
                toc;
            end
        end
    end
end

corrFile = sprintf('%s_corr%s',figPrefix,figSuffix);
save(corrFile,'R','rho','prho','chi2','pchi2');

cv = sqrt(cv2);

%%

indices = triu(true(n),1);

dataMs = {rho.^2 prho chi2 pchi2 cv cv2};
dataVs = zeros(N,8,4,6);
prctiles = zeros(7,8,4,6);

for hh = 1:6
    figure;
    
    for ii = 1:4
        subplot(2,2,ii);

        dataV = zeros(N,8);

        for jj = 1:8
            dataM = datas{hh}(:,:,jj,ii);
            dataV(:,jj) = dataM(indices);
        end

        boxplot(dataV);
        
        prctiles(:,:,ii,hh) = prctile(dataV,[1 5 25 50 75 95 99]);
        
        dataVs(:,:,ii,hh) = dataV;
    end
end

%%

disp(squeeze(prctiles(6,:,:,1)));
disp(squeeze(prctiles(2,:,:,2)));
disp(squeeze(prctiles(6,:,:,3)));
disp(squeeze(prctiles(2,:,:,4)));
disp(squeeze(prctiles(6,:,:,5)));
disp(squeeze(prctiles(6,:,:,6)));

%%

disp(sum(reshape(dataVs(:,:,:,2),N*32,1) < 0.05)/(N*32));
disp(sum(reshape(dataVs(:,:,:,2),N*32,1) < 0.01)/(N*32));
disp(sum(reshape(dataVs(:,:,:,2),N*32,1) < 0.001)/(N*32));
disp(sum(reshape(dataVs(:,:,:,2),N*32,1) < 0.05/(N*32))/(N*32));
disp(sum(fdrcorrect(reshape(dataVs(:,:,:,2),N*32,1)))/(N*32));
disp(sum(holmBonferroni(reshape(dataVs(:,:,:,2),N*32,1)))/(N*32));

%%

disp(sum(reshape(dataVs(:,:,:,4),N*32,1) < 0.05)/(N*32));
disp(sum(reshape(dataVs(:,:,:,4),N*32,1) < 0.01)/(N*32));
disp(sum(reshape(dataVs(:,:,:,4),N*32,1) < 0.001)/(N*32));
disp(sum(reshape(dataVs(:,:,:,4),N*32,1) < 0.05/(N*32))/(N*32));
disp(sum(fdrcorrect(reshape(dataVs(:,:,:,4),N*32,1)))/(N*32));
disp(sum(holmBonferroni(reshape(dataVs(:,:,:,4),N*32,1)))/(N*32));

save(corrFile,'dataMs','dataVs','prctiles','-append');

%%

assert(false);

%% COUNT

ws = warning('off','GetOpt:UnknownOption');
pids = biasCorrectedPID(stimuli(trials,:),[resp1 resp2],'savefile',[figPrefix '_synergy_temp']);
warning(ws);
% pids = zeros(N,4,4);
% 
% tic;
% 
% for hh = 1:4
%     oo = 1;
%     for ii = 1:n-1
%         for jj = (ii+1):n
%             nx1 = max(responses{indices(ii)}(:,1,1,hh))+1;
%             nx2 = max(responses{indices(jj)}(:,1,1,hh))+1;
%             counts = zeros(8,nx1,nx2);
%             for kk = 1:size(stimuli,1)
%                 ll = stimuli(kk,1);
%                 mm = responses{indices(ii)}(kk,1,1,hh)+1;
%                 nn = responses{indices(jj)}(kk,1,1,hh)+1;
%                 counts(ll,mm,nn) = counts(ll,mm,nn) + 1;
%             end
%             pids(oo,:,hh) = PID(counts,sources,sets,trans)';
%             oo = oo + 1;
%         end
%     end
% end
% 
% pids(:,5,:) = sum(pids,2);
% 
% toc;

% tic;
% 
% % TODO : test. Also pretty sure that's not a valid indexing expression.
% pids = pairwisePID(stimuli,{responses{:}(:,1,1,:)});
% 
% toc;

%%

if isMCS
    figInfix = '';
else
    figInfix = '25';
end

filename = sprintf('%s_synergy%s%s',figPrefix,figInfix,figSuffix);
save(filename,'pids');

%%
pids(any(any(abs(pids) > 1e10,3),2),:,:) = NaN;

figure;
barwitherr(squeeze(nanstd(pids))/sqrt(N),squeeze(nanmean(pids)));
ylabel('Information (bits)')
set(gca,'XTickLabel',{'Redundancy' 'Unique 1' 'Unique 2' 'Synergy' 'Information'})
legend(legends,'Location','NorthWest')

%%

saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

%%

% save(filename,'pids','sources','sets','trans');

%% LATENCY
% 
% bws = [100 50 25 10 5 2.5 1];
% nbs = 500./bws + 1;
% multipliers = 1000./bws;
% 
% %%
% 
% % pidl = cell(7,1);
% 
% for gg = 6:7
%     pidl{gg} = zeros(N,4,4);
% 
%     tic;
%     
%     for hh = 1:4
%         oo = 1;
%         for ii = 1:n-1
%             for jj = (ii+1):n
% %                 tic;
%                 counts = zeros(8,nbs(gg),nbs(gg));
%                 for kk = 1:size(stimuli,1)
%                     nn = stimuli(kk,1);
% 
%                     resps = [responses{indices(ii)}(kk,2,1,hh) responses{indices(jj)}(kk,2,1,hh)];
%                     mm = [0 0];
% 
%                     for ll = 1:2
%                         resp = resps(ll);
% 
%                         if resp > 0.5
%                             mm(ll) = nbs(gg);
%                         elseif resp == 0.5
%                             mm(ll) = nbs(gg)-1;
%                         else
%                             mm(ll) = floor(resp*multipliers(gg))+1;
%                         end
%                     end
% 
%                     counts(nn,mm(1),mm(2)) = counts(nn,mm(1),mm(2)) + 1;
%                 end
%     %             tic;
%                 pidl{gg}(oo,:,hh) = PID(counts,sources,sets,trans)';
%     %             toc;
%                 oo = oo + 1;
% %                 toc;
%             end
%         end
%     end
% 
%     pidl{gg}(:,5,:) = sum(pidl{gg},2);
% 
%     toc;
% 
%     %%
% 
%     barwitherr(squeeze(std(pidl{gg}))/sqrt(N),squeeze(mean(pidl{gg})));
%     ylabel('Information (bits)')
%     set(gca,'XTickLabel',{'Redundancy' 'Unique 1' 'Unique 2' 'Synergy' 'Information'})
%     legend(legends,'Location','NorthWest')
% 
%     %%
%     figFile = sprintf('%s_synergy25_latency_%d_bins',figPrefix,nbs(gg));
%     saveas(gcf,figFile,'fig');
%     saveas(gcf,figFile,'png');
%     
%     save('synergy25.mat','pidl','-append')
% end

%% SPIKE TRAIN

% if isMCS
%     figInfix = '';
% else
%     figInfix = '25';
% end
% 
% filename = sprintf('%s_synergy%s%s',figPrefix,figInfix,figSuffix);
% 
% if exist(filename,'file')
%     load(filename,'pidt','info1','info2','kld1');
%     
%     pairs = zeros(N,2);
%     
%     nn = 0;
%     for ii = 1:n-1
%         for jj = ii+1:n
%             nn = nn + 1;
% 
%             pairs(nn,:) = [ii jj];
%         end
%     end
% else
%     load('sip.mat', 'stimuli', 'responses', 'responseIndices')
% 
%     %%
% 
%     info1 = zeros(n,4);
%     kld1 = zeros(n,8,4);
% 
%     for ii = 1:n
%         for jj = 1:4
%             tic;
%             ts = squeeze(responses{indices(ii)}(1,trials,1,jj))';
%             ss = squeeze(responses{indices(ii)}(2,trials,1,jj))';
%             es = squeeze(responses{indices(ii)}(3,trials,1,jj))';
%             r = arrayfun(@(t,s,e) spiketimestamps{idx(ii)-1}(s:e)-t,ts,ss,es,'UniformOutput',false);
%             kld1(ii,:,jj) = binlessInfo(stimuli(trials,jj),r,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1,'kld',true); %#ok<NODEF>
%             info1(ii,jj) = binlessInfo(stimuli(trials,jj),r,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1);
%             toc;
%         end
%     end
% 
%     info1 = max(info1,0);
% 
%     info2 = zeros(N,4);
%     pairs = zeros(N,2);
%     nTrials = numel(trials);
% 
%     nn = 0;
%     for ii = 1:n-1
%         for jj = ii+1:n
%             nn = nn + 1;
% 
%             pairs(nn,:) = [ii jj];
%             T = spiketimestamps(idx([ii jj])-1);
%             R = cat(5,responses{indices([ii jj])});
% 
%             for kk = 1:4
%                 tic;
%                 ts = squeeze(R(1,trials,1,kk,:));
%                 ss = squeeze(R(2,trials,1,kk,:));
%                 es = squeeze(R(3,trials,1,kk,:));
%                 r = arrayfun(@(t,s,e,c) T{c}(s:e)-t,ts,ss,es,[ones(nTrials,1) 2*ones(nTrials,1)],'UniformOutput',false);
%                 info2(nn,kk) = binlessInfo(stimuli(trials,kk),r,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1);
%                 toc;
%             end
%         end
%     end
% 
%     info2 = max(info2,0);
% end
% 
% %%
% 
% red = squeeze(mean(min(cat(4,kld1(pairs(:,1),:,:),kld1(pairs(:,2),:,:)),[],4),2));
% uq1 = info1(pairs(:,1),:,1)-red;
% uq2 = info1(pairs(:,2),:,1)-red;
% syn = info2-uq1-uq2-red;
% pidt = permute(cat(3,red,uq1,uq2,syn,info2),[1 3 2]);
% 
% save(filename,'pidt','info1','info2','kld1');
% 
% %%
% pids(any(any(abs(pids) > 1e10,3),2),:,:) = NaN;
% pidt(any(any(abs(pidt) > 1e10,3),2),:,:) = NaN;
% %%
% figure;
% % barwitherr(squeeze(nanstd(pidt))/sqrt(N),squeeze(nanmean(pidt)));
% barwitherr(cat(3,squeeze(prctile(pidt,75))-squeeze(nanmedian(pidt)),squeeze(nanmedian(pidt))-squeeze(prctile(pidt,25))),squeeze(nanmedian(pidt)));
% ylabel('Information (bits)')
% set(gca,'XTickLabel',{'Redundancy' 'Unique 1' 'Unique 2' 'Synergy' 'Information'})
% legend(legends,'Location','NorthWest')
% 
% %%
% 
% saveas(gcf,filename,'fig');
% saveas(gcf,filename,'png');