if false
expts = {   ...
    'P95_P96Thy1\P96',  ...
    47,                 ...
    48,                 ...
    49,                 ...
    50,                 ...
    54,                 ...
    66,                 ...
    69                  ...
    };

stims = {       ...
    [1 2 3 5],  ...
    [1 2 3 5 6],    ...
    [1 2 3 5 6],    ...
    [1 3 3 5 6],    ...
    [1 2 3 5 6],    ...
    [1 2 3 5 6],    ...
    [1 2 4 5 8 9],    ...
    [2 3 4 5 6 7 8 10]   ...
    };

stimulusFiles = { ...
    NaN ...
    {NaN 'contrast sequence' 'frequency sequence' 'contrast sequence' 'frequency sequence'} ... 
    {NaN 'contrast sequence' 'frequency sequence' 'contrast sequence' 'frequency sequence'} ... 
    {NaN 'contrast sequence' 'frequency sequence' 'contrast sequence' 'frequency sequence'} ... 
    {NaN 'contrast sequence' 'frequency sequence' 'contrast sequence' 'frequency sequence'} ... 
    {NaN 'contrast sequence' 'frequency sequence' 'contrast sequence' 'frequency sequence'} ... 
    {'frequency sequence' 'lumstep sequence' 'frequency sequence' 'lumstep sequence' 'lumstep sequence' 'frequency sequence'} ...
    {'lumstep sequence' 'lumstep sequence' 'constep sequence' 'revgrat sequence' 'lumstep sequence' 'constep sequence' 'revgrat sequence' 'lumstep sequence'} ...
    };

detectionFiles = {  ...
    {'ChR2_cells_ff' 'ChR2_cells_ff' 'ChR2_cells_ff' 'ChR2_cells_ff'} ...
    {'ChR2_cells_ff' 'ChR2_cells_contrast_ctrl' 'ChR2_cells_frequency_ctrl' 'ChR2_cells_contrast_drug' 'ChR2_cells_frequency_drug'}  ...
    {'ChR2_cells_ff' 'ChR2_cells_contrast_ctrl' 'ChR2_cells_frequency_ctrl' 'ChR2_cells_contrast_drug' 'ChR2_cells_frequency_drug'}  ...
    {'ChR2_cells_ff' 'ChR2_cells_contrast_ctrl' 'ChR2_cells_frequency_ctrl' 'ChR2_cells_contrast_drug' 'ChR2_cells_frequency_drug'}  ...
    {'ChR2_cells_ff' 'ChR2_cells_contrast_ctrl' 'ChR2_cells_frequency_ctrl' 'ChR2_cells_contrast_drug' 'ChR2_cells_frequency_drug'}  ...
    {'ChR2_cells_ff' 'ChR2_cells_contrast_ctrl' 'ChR2_cells_frequency_ctrl' 'ChR2_cells_contrast_drug' 'ChR2_cells_frequency_drug'}  ...
    {'ChR2_cells_frequency_ctrl' 'ChR2_cells_lumstep_ctrl' 'ChR2_cells_frequency_drug' 'ChR2_cells_lumstep_drug' 'ChR2_cells_lumstep_high' 'ChR2_cells_frequency_high'}  ...
    {'ChR2_cells_lumstep_lowc' 'ChR2_cells_lumstep_ctrl' 'ChR2_cells_constep_ctrl' 'ChR2_cells_revgrat_ctrl' 'ChR2_cells_lumstep_drug' 'ChR2_cells_constep_drug' 'ChR2_cells_revgrat_drug' 'ChR2_cells_lumstep_lowd'}  ...
    };

spontStartIndices = {   ...
    [16 16 16 16] ...
    [16 1 1 1 1] ...
    [16 1 1 1 1] ...
    [16 1 1 1 1] ...
    [16 1 1 1 1] ...
    [16 1 1 1 1] ...
    [1 1 1 1 1 1] ...
    [1 1 1 1 1 1 1 1] ...
    };

whiteConditions = {
    NaN ...
    [NaN 49 49 49 49] ...
    [NaN 49 49 49 49] ...
    [NaN 49 49 49 49] ...
    [NaN 49 49 49 49] ...
    [NaN 49 49 49 49] ...
    [49 9 49 9 9 49] ...
    [9 9 26 7 9 26 7 9] ...
    };

%%

for ii = 6:8
    if ischar(expts{ii})
        cd(expts{ii})
    else
        cd(sprintf('JBOG%04d',expts{ii}));
    end
    
    if ii == 1
        detectChR2ResponsiveCells(stims{ii},0,0,0,5,'','',false,false,false,30);
    else
        for jj = (1+4*(ii==6)):numel(stims{ii})
            if strcmp(detectionFiles{ii}{jj},'ChR2_cells_ff')
                detectChR2ResponsiveCells(stims{ii}(jj),0,0,0,5,'','',false,false,false,30);
            elseif ii == 6 && jj == 5
                detectChR2ResponsiveCells(stims{ii}(jj),0,whiteConditions{ii}(jj),4,5,[stimulusFiles{ii}{jj} '.mat'],[detectionFiles{ii}{jj}(11:end) 'new'],true,false,false,47);
            else
                detectChR2ResponsiveCells(stims{ii}(jj),0,whiteConditions{ii}(jj),4,5,[stimulusFiles{ii}{jj} '.mat'],[detectionFiles{ii}{jj}(11:end) 'new'],true,false,false);
            end
        end
    end
    
    cd ..;
    
    if ii == 1
        cd ..
    end
end

%%

for hh = 1:8
    if hh == 1
        cd(expts{hh});
    else
        cd(sprintf('JBOG%04d',expts{hh}));
    end
    
    for gg = 1:numel(stims{hh})
        detectionFile = detectionFiles{hh}{gg};
        
        if strcmp(detectionFile,'ChR2_cells_ff')
            saveFile = 'ChR2_cells_new';
        else
            detectionFile = ['ChR2_cells_' detectionFile(11:end) 'new'];
            saveFile = detectionFile;
            saveFile(11) = [];
            saveFile(end-2:end+1) = '_new';
        end
        
        spontStartIndex = spontStartIndices{hh}(gg);
        
        newDetection;
    end
    
    cd ..
    
    if hh == 1
        cd ..
    end
end

%%

kcls = {        ...
    [3 3 6 9],  ...
    [3 9 9 9 9],    ...
    [3 9 9 3 3],    ...
    [3 9 9 9 9],    ...
    [3 9 9 9 9],    ...
    [3 9 9 9 9],    ...
    [3 3 3 3 9 9],    ...
    [3 9 9 9 9 9 9 3]   ...
    };

mfas = {                        ...
    [false true(1,3)],          ...
    [false false false true true],         ...
    [false false false true true],         ...
    [false false false true true],         ...
    [false false false true true],         ...
    [false false false true true],         ...
    [false false true true true true],          ...
    [false false false false true true true true]     ...
    };
end

%%

for hh = 6:8
    if hh == 1
        P = 96;
        cd(expts{hh});
    else
        P = sprintf('JBOG%04d',expts{hh});
        cd(P);
    end
    
    detectionFile = detectionFiles{hh};
    
    stimIndices = cell(numel(detectionFile),2);
    
    for gg = 1:numel(detectionFile)
        if strcmp(detectionFile{gg}(end-1:end),'ff')
            detectionFile{gg} = [detectionFile{gg}(1:end-2) 'new.mat'];
            stimIndices{gg,1} = 1:2:59;
            stimIndices{gg,2} = 2:2:60;
        else
            detectionFile{gg} = [detectionFile{gg} '_new.mat'];
            load(stimulusFiles{hh}{gg});
            stimIndices{gg,1} = 4*(find(conditionOrder(:,1) == whiteConditions{hh}(gg))-1)+1;
            stimIndices{gg,2} = 4*(find(conditionOrder(:,1) == 0)-1)+1;
        end
    end
    
    if hh == 6
        stimIndices{5,1} = stimIndices{5,1}(1:47);
        stimIndices{5,2} = stimIndices{5,2}(1:47);
    end
        
    recordings = stims{hh};
    
    flab = 'MFA';
    xlab = 'KCl Concentration (mM)';
    xticks = kcls{hh};
    xfill = [repmat(find(mfas{hh})-0.5,2,1); repmat(find(mfas{hh})+0.5,2,1)];
    yfill = repmat([0.1 40 40 0.1]',1,size(xfill,2));
    
    rerereanalyseBasicChR2Experiments;
    
    cd ..
    
    if hh == 1
        cd ..
    end
end

%%

evokedFRs = cell(size(expts));
nCells = cell(size(expts));

for ii = 1:8
    expt = expts{ii};
    
    if ischar(expt)
        load(sprintf('%s\\P96_responses_new.mat',expt));
        cd(expt)
    else
        load(sprintf('JBOG%04d\\JBOG%04d_responses_new.mat',expt,expt));
        cd(sprintf('JBOG%04d',expt));
    end
    
    evokedFRs{ii} = mdfr;
    
    detectionFile = detectionFiles{ii};
    
    load([strrep(detectionFile{1},'_ff','') '_new.mat'],'allBestUnits');
    nCells{ii} = cellfun(@numel,allBestUnits);
    
    if numel(unique(detectionFile)) == 1
        cd ..
        cd ..
        continue;
    end
    
    for jj = 2:numel(detectionFile)
        load([strrep(detectionFile{jj},'_ff','') '_new.mat'],'allBestUnits');
        nCells{ii}(jj,:) = cellfun(@numel,allBestUnits);
    end
    
    cd ..
end

%%

E = vertcat(evokedFRs{:});
E = E(:);

N = vertcat(nCells{:});
N = N(:);

%%

topDir = 'D:\John B\Electrophysiology\Optogenetics\';

cd(topDir);

cleanup = onCleanup(@() cd(topDir));

%%

nExperiments = numel(expts);
spontTimes = cell(nExperiments,1);
nSpikess = cell(nExperiments,1);
sampleRate = 7055.258405732180;

for ii = 1:nExperiments
    expt = expts{ii};
    
    if ischar(expt)
        cd(expt);
    else
        cd(sprintf('JBOG%04d',expt));
    end
    
%     nStims = numel(stims{ii});
%     nSamples = zeros(nStims,1);
    fin = fopen('originalSplit.txt');
    closeFile = onCleanup(@() fclose(fin));
    
%     jj = 0;
%     while ~feof(fin)
%         data = textscan(fin, '%*s %d', 'CommentStyle', '#', 'Delimiter', ',');
%         
%         if isempty(data)
%             continue
%         end
%         
%         jj = jj + 1;
%         nSamples(jj) = data{jj};
%     end

    nSamples = textscan(fin, '%*s %d', 'CommentStyle', '#', 'Delimiter', ',');
    nSamples = nSamples{1};
    
    stim = stims{ii};
    nStims = numel(stim);
    assert(numel(nSamples) >= nStims);
    
    t = load('EventNo.mat','EventNo');
    startTimes = cellfun(@(t) t(1),t.EventNo(stim))';
    endTimes = [0; cumsum(double(nSamples)/sampleRate)];
    endTimes = endTimes(stim);
    
    assert(all(startTimes > endTimes));
    
    c = load('channelNames.mat','channelNames');
    s = load('spiketimestamps.mat','spiketimestamps');
    spikeTimes = s.spiketimestamps([c.channelNames{6,:}] == 0);
    n = numel(spikeTimes);
    nSpikess{ii} = zeros(n,2,2);
    spontTimes{ii} = zeros(1,2,2);
    
    for jj = 1:nStims
        if kcls{ii}(jj) == 6
            continue
        end
        
        spontTime = startTimes(jj)-endTimes(jj);
        
        kidx = (kcls{ii}(jj)+3)/6;
        midx = mfas{ii}(jj)+1;
        
        spontTimes{ii}(1,kidx,midx) = spontTimes{ii}(1,kidx,midx)+spontTime;
        
        nSpikess{ii}(:,kidx,midx) = nSpikess{ii}(:,kidx,midx)+cellfun(@(s) sum(s > endTimes(jj) & s <= startTimes(jj)),spikeTimes)';
    end
    
    cd(topDir);
end

%%

firingRates = cellfun(@(f,t) f./repmat(t,size(f,1),1),nSpikess,spontTimes,'UniformOutput',false);

%%

K = [kcls{:}];
M = [mfas{:}];

%%

figure;
boxplot(E(repmat(K ~= 6,2,1)),{repmat(40*M(K ~= 6)',2,1) repmat(K(K ~= 6)',2,1) kron([0;1],ones(sum(K ~= 6),1))})

figure;
boxplot(N(repmat(K ~= 6,2,1)),{repmat(40*M(K ~= 6)',2,1) repmat(K(K ~= 6)',2,1) kron([0;1],ones(sum(K ~= 6),1))})

%%

FR = cellfun(@median,firingRates,'UniformOutput',false);
FR = [FR{:}];

%%

uK = unique(K);
nK = numel(uK);

uM = unique(M);
nM = numel(uM);

nFR = zeros(nK,nM);
mFR = nan(nK,nM);
sFR = nan(nK,nM);

for ii = 1:nK
    for jj = 1:nM
        fr = FR(K == uK(ii) & M == uM(jj));
        disp(fr);
        nFR(ii,jj) = numel(fr);
        mFR(ii,jj) = mean(fr);
        sFR(ii,jj) = std(fr);
    end
end

%%

figure;
boxplot(FR(K ~= 6),{40*M(K ~= 6) K(K ~= 6)});
set(gca,'XTickLabel',{' '});
ylabel('Spontaneous Firing Rate (Hz)');
yy = ylim;
ylim([-0.1 4.5]);

mfaLabels = arrayfun(@(c) sprintf('%d {\\mu}M MFA\n',c), [0 40], 'UniformOutput', false);
kclLabels = arrayfun(@(c) sprintf('+ %d mM KCl', c), [3 9], 'UniformOutput', false);

repLabelHeights = [0.5 0.1; 0.4 0.95];
repLabelAlignments = {'middle' 'bottom'; 'cap' 'top'};

for ii = 1:2
    for jj = 1:2
        label = [mfaLabels{ii} kclLabels{jj}];
        x = jj+2*(ii-1);
        text(x,-0.2,label,'HorizontalAlignment','center','VerticalAlignment','top');
        text(x,repLabelHeights(ii,jj),sprintf('n = %d',sum(M == uM(ii) & K == uK(2*jj-1))),'HorizontalAlignment','center','VerticalAlignment',repLabelAlignments{ii,jj});
    end
end

%%

disp(ranksum(FR(M == uM(1) & K == uK(1)),FR(M == uM(1) & K == uK(3))))
disp(ranksum(FR(M == uM(2) & K == uK(1)),FR(M == uM(2) & K == uK(3))))

%%

save('kclvspont.mat','kcls','mfas','firingRates','evokedFRs','nCells','K','M','FR','E','N');