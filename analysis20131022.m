%%
dirs = { ...
    'JBWT0024', ...
    'JBWT0025', ...
    'JBWT0026', ...
    'JBWT0029', ...
    'JBWT0032', ...
    'JBWT0033', ...
    'JBWT0034', ...
    'JBWT0035', ...
    'JBWT0036' ...
    };
sigs = {'visual', 'visual', 'visual', 'rf', 'visual', 'visual', 'visual', 'visual', 'visual', 'visual'};
recs = {[2 4 7 9], [2 4], 4, [2 4], [4 5], [4 5], [4 5], [3 4], [3 4]};
ages = [26 99 12 73 19 23 31 26 46];
lights = {[0 1 0 1],[0 1],1,[0 1],[0 0],[0 0],[0 0],[0 0],[0 0]};
sides = {'llrr','rr','u','rr','ll','ll','uu','ll','ll'};
tasks = {'ffff','ff','f','ff','fc','fc','fc','fc','fc'};
% dirs = { ...
%     'JBWT0024', ...
%     'JBWT0025', ...
%     'JBWT0029', ...
%     'JBWT0032', ...
%     'JBWT0033', ...
%     'JBWT0034', ...
%     'JBWT0035', ...
%     'JBWT0036' ...
%     };
% sigs = {'visual', 'visual', 'rf', 'visual', 'visual', 'visual', 'visual', 'visual', 'visual'};
% recs = {[2 4 7 9], [2 4], [2 4], [4 5], [4 5], 5, [3 4], [3 4]};
% ages = [26 99 73 19 23 31 26 46];
% lights = {[0 1 0 1],[0 1],[0 1],[0 0],[0 0],0,[0 0],[0 0]};
% sides = {'llrr','rr','rr','ll','ll','u','ll','ll'};
% tasks = {'ffff','ff','ff','fc','fc','c','fc','fc'};
% dirs = { ...
%     'JBWT0024', ...
%     'JBWT0025', ...
%     'JBWT0026', ...
%     'JBWT0029', ...
%     'JBWT0032', ...
%     'JBWT0034', ...
%     'JBWT0036' ...
%     };
% sigs = {'visual', 'visual', 'visual', 'rf', 'visual', 'visual', 'visual', 'visual', 'visual'};
% recs = {[2 4 7 9], [2 4], 4, [2 4], [4 5], [4 5], [3 4]};
% ages = [26 99 12 73 19 31 46];
% lights = {[0 1 0 1],[0 1],1,[0 1],[0 0],[0 0],[0 0]};
% sides = {'llrr','rr','u','rr','ll','uu','ll'};
% tasks = {'ffff','ff','f','ff','fc','fc','fc'};
%%
perf = [];
nneurons = [];
all = [];
code = [];
width = [];
contrast = [];
polarity = [];
age = [];
light = [];
side = [];
task = [];
retina = [];
superTotalNeurons = 0;

currentDir = pwd;

r = 1;
for ii = 1:numel(dirs)
    try
        cd(dirs{ii});
        rec = recs{ii};
        
        load('./initRecordings.mat');
        
        for jj = 1:numel(rec)
            if jj == 3
                r = r + 1;
            end
            
            recording = recordings(rec(jj));
            fileDir = getAnalysisOutputDir(recording);
            load(sprintf('%s\\%s_bayes_bars_1_peak_50',fileDir,recording.dataFile));
            
            perfs = zeros(0,4,2,4);
            ns = zeros(0,4,2,4);
            as = zeros(0,4,2,4);
            
            for kk = 1:numel(performances)
                performance = performances{kk};
                
                if kk == 1
                    totalNeurons = size(performance,1);
                    superTotalNeurons = superTotalNeurons + totalNeurons;
                    nNeurons = 1;
                elseif kk == numel(performances)
                    nNeurons = totalNeurons;
                else
                    nNeurons = 5*(kk-1);
                end
                
                perfs(end+1:end+size(performance,1),:,:,:) = performance;
                ns(end+1:end+size(performance,1),:,:,:) = nNeurons;
                as(end+1:end+size(performance,1),:,:,:) = kk == numel(performances);
            end
               
            n = size(perfs,1);
            
            p = reshape(perfs,n*4*2*4,1);
            perf = [perf; p]; %#ok<AGROW>
            nneurons = [nneurons; reshape(ns,n*4*2*4,1)]; %#ok<AGROW>
            all = [all; reshape(as,n*4*2*4,1)]; %#ok<AGROW>
            
            c = reshape(1:4,[1 1 1 4]);
            c = repmat(c,[n 4 2 1]);
            c = reshape(c,n*4*2*4,1);
            code = [code; c]; %#ok<AGROW>
            
            if tasks{ii}(jj) == 'f'
                % this is wrong for one recording but eh
                w = reshape(repmat(2.^(3:6),[n 1 2 4]),size(c));
                k = ones(size(c));
            else
                w = 40*ones(size(c));
                k = reshape(repmat([1.0 0.4 0.2 0.1],[n 1 2 4]),size(c));
            end
            
            p = zeros(1,1,2);
            p(1,1,1) = 1;
            p(1,1,2) = 2;
            
            p = reshape(repmat(p,[n 4 1 4]),size(c));
            
            width = [width; w]; %#ok<AGROW>
            contrast = [contrast; k]; %#ok<AGROW>
            polarity = [polarity; p]; %#ok<AGROW>
            
            age = [age; ages(ii)*ones(size(p))]; %#ok<AGROW>
            light = [light; lights{ii}(jj)*ones(size(p))]; %#ok<AGROW>
            side = [side; repmat(sides{ii}(jj),size(p))]; %#ok<AGROW>
            task = [task; repmat(tasks{ii}(jj),size(p))]; %#ok<AGROW>
            retina = [retina; r*ones(size(p))]; %#ok<AGROW>
        end
        
        r = r + 1;
        cd ..;
    catch err
        logMatlabError(err);
        cd(currentDir);
    end
end


%%
save('aggregate_bayes','perf','nneurons','code','width','contrast','polarity','age','light','side','task','retina');

%%
ns = unique(nneurons);
ns = [1; ns(ns < 40 & ceil(ns/5) == floor(ns/5))];

mperf = zeros(numel(ns),2,4);
lperf = zeros(numel(ns),2,4);
uperf = zeros(numel(ns),2,4);

coeffs = zeros(3,2,4);
% coeffs = zeros(2,2,4);
resids = cell(2,2,4);
sigmas = cell(2,4);

for ii = 1:2
    for jj = 1:4
        for kk = 1:numel(ns)
            index = polarity == ii & code == jj & nneurons == ns(kk);
            p = perf(index);
            
            if isempty(p)
                continue;
            end
            
            mperf(kk,ii,jj) = median(p);
            lperf(kk,ii,jj) = prctile(p,25);
            uperf(kk,ii,jj) = prctile(p,75);
        end
        
        index = polarity == ii & code == jj;
        y = perf(index);
        x = nneurons(index);
        
%         for kk = 1:10
%             [q,S] = polyfit(x,y,kk);
%             polyfits{kk,ii,jj} = q;
%             polystats{kk,ii,jj} = S;
%         end

        [b,R,~,S] = nlinfit(x,y,@(b,x) 1*(1-exp(-b*x)),1);
        coeffs(1,ii,jj) = b;
        resids{1,ii,jj} = R;
        sigmas{ii,jj} = S;
        
        [b,R,~,S] = nlinfit(x,y,@(b,x) b(1)*(1-exp(-b(2)*x)),[1 1]);
        coeffs(2:3,ii,jj) = b;
        resids{2,ii,jj} = R;
        sigmas{ii,jj} = S;
    end
end

lperf = mperf-lperf;
uperf = uperf-mperf;

%% 

cis = zeros(2,2,4);

for ii = 1:2
    for jj = 1:4
        cis(:,ii,jj) = nlparci(coeffs(ii,jj),resids{ii,jj},'covar',sigmas{ii,jj});
    end
end

%%

figure;
colours = 'rgbm';
codes = {'Count' 'Latency' 'Timing' 'Correlation'};

% index = [1; find(ns < 35 & ceil(ns/5) == floor(ns/5))];
x = 0:5:max(ns);

for ii = 1
%     subplot(1,2,ii);
    hold on;
    
    for jj = 1:3
%         errorbar(ns(index),mperf(index,ii,jj),lperf(index,ii,jj),uperf(index,ii,jj),'Color',colours(jj));
        f = @(x) 1-exp(-coeffs(1,ii,jj)*x);
%         f = @(x) coeffs(2,ii,jj)*(1-exp(-coeffs(3,ii,jj)*x));
        errorbar(ns,mperf(:,ii,jj),lperf(:,ii,jj),uperf(:,ii,jj),'Color',colours(jj));
        [X,Y] = fplot(f,xlim); 
        plot(X,Y,'Color',colours(jj),'LineStyle','--');
    end
end

% legend(codes,'Location','SouthEast');

%%

widths = unique(width);
widths = widths([1 2 3 5]);
contrasts = unique(contrast);

mpvw = zeros(numel(widths),2,4);
lpvw = zeros(numel(widths),2,4);
upvw = zeros(numel(widths),2,4);

mpvc = zeros(numel(contrasts),2,4);
lpvc = zeros(numel(contrasts),2,4);
upvc = zeros(numel(contrasts),2,4);

for ii = 1:2
    for jj = 1:4
        for kk = 1:numel(widths)
            index = polarity == ii & code == jj & width == widths(kk) & nneurons == 1;
            p = perf(index);
            
            if isempty(index)
                continue;
            end
            
            mpvw(kk,ii,jj) = median(p);
            lpvw(kk,ii,jj) = prctile(p,25);
            upvw(kk,ii,jj) = prctile(p,75);
        end
        
        for kk = 1:numel(contrasts)
            index = polarity == ii & code == jj & contrast == contrasts(kk) & nneurons == 1;
            p = perf(index);
            
            if isempty(index)
                continue;
            end
            
            mpvc(kk,ii,jj) = median(p);
            lpvc(kk,ii,jj) = prctile(p,25);
            upvc(kk,ii,jj) = prctile(p,75);
        end
    end
end

lpvw = mpvw-lpvw;
upvw = upvw-mpvw;

lpvc = mpvc-lpvc;
upvc = upvc-mpvc;

%%

lineStyles = {'-' '--'};
figure;
hold on;

for ii = 1:2
    for jj = 1:3
        errorbar(widths,mpvw(:,ii,jj),lpvw(:,ii,jj),upvw(:,ii,jj),'Color',colours(jj),'LineStyle',lineStyles{ii});
    end
end

%%

figure;
hold on;

for ii = 1:2
    for jj = 1:3
        errorbar(contrasts,mpvc(:,ii,jj),lpvc(:,ii,jj),upvc(:,ii,jj),'Color',colours(jj),'LineStyle',lineStyles{ii});
    end
end

%%

mapvw = zeros(numel(widths),2,4);
lapvw = zeros(numel(widths),2,4);
uapvw = zeros(numel(widths),2,4);

mapvc = zeros(numel(contrasts),2,4);
lapvc = zeros(numel(contrasts),2,4);
uapvc = zeros(numel(contrasts),2,4);

for ii = 1:2
    for jj = 1:4
        for kk = 1:numel(widths)
            index = polarity == ii & code == jj & width == widths(kk) & all == 1;
            p = perf(index);
            
            if isempty(index)
                continue;
            end
            
            mapvw(kk,ii,jj) = median(p);
            lapvw(kk,ii,jj) = prctile(p,25);
            uapvw(kk,ii,jj) = prctile(p,75);
        end
        
        for kk = 1:numel(contrasts)
            index = polarity == ii & code == jj & contrast == contrasts(kk) & all == 1;
            p = perf(index);
            
            if isempty(index)
                continue;
            end
            
            mapvc(kk,ii,jj) = median(p);
            lapvc(kk,ii,jj) = prctile(p,25);
            uapvc(kk,ii,jj) = prctile(p,75);
        end
    end
end

lapvw = mpvw-lapvw;
uapvw = upvw-mapvw;

lapvc = mapvc-lapvc;
uapvc = uapvc-mapvc;

%%

lineStyles = {'-' '--'};
figure;
hold on;

for ii = 1:2
    for jj = 1:3
        errorbar(widths,mapvw(:,ii,jj),lapvw(:,ii,jj),uapvw(:,ii,jj),'Color',colours(jj),'LineStyle',lineStyles{ii});
    end
end

%%

figure;
hold on;

for ii = 1:2
    for jj = 1:3
        errorbar(contrasts,mapvc(:,ii,jj),lapvc(:,ii,jj),uapvc(:,ii,jj),'Color',colours(jj),'LineStyle',lineStyles{ii});
    end
end

%%

sperfw = zeros(max(retina),4,numel(widths),2);
aperfw = zeros(max(retina),4,numel(widths),2);
sperfc = zeros(max(retina),4,numel(contrasts),2);
aperfc = zeros(max(retina),4,numel(contrasts),2);

for ii = 1:max(retina)
    for jj = 1:4
        for kk = 1:2
            index1 = retina == ii & code == jj & polarity == kk;
            
            for ll = 1:numel(widths)
                index2 = index1 & task == 'f' & width == widths(ll);
                sperfw(ii,jj,ll,kk) = median(perf(index2 & nneurons == 1));
                aperfw(ii,jj,ll,kk) = perf(index2 & all == 1 & light == (ii == 4));
            end
        
            for ll = 1:numel(contrasts)
                index2 = index1 & task == 'c' & contrast == contrasts(ll);
                
                if sum(index2) == 0
                    continue;
                end
                
                sperfc(ii,jj,ll,kk) = median(perf(index2 & nneurons == 1));
                aperfc(ii,jj,ll,kk) = perf(index2 & all == 1);
            end
        end
    end
end

sperfc = sperfc(6:end,:,:,:);
aperfc = aperfc(6:end,:,:,:);

%%

for ii = 1:2
    data = reshape(sperfw(:,:,:,ii),size(sperfw,1)*4,numel(widths));
    friedman(data,size(data,1)/4);
    
    data = reshape(aperfw(:,:,:,ii),size(aperfw,1)*4,numel(widths));
    friedman(data,size(data,1)/4);
    
    data = reshape(sperfc(:,:,:,ii),size(sperfc,1)*4,numel(contrasts));
    friedman(data,size(data,1)/4);
    
    data = reshape(aperfc(:,:,:,ii),size(aperfc,1)*4,numel(contrasts));
    friedman(data,size(data,1)/4);
end