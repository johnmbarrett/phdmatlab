dirs = { ...
    'JBWT0032', ...
    'JBWT0033', ...
    'JBWT0034', ...
    'JBWT0035', ...
    'JBWT0036' ...
%     'F:\John B\Electrophys\Wild-Type\JBWT0026' ...
%     '../JBWT0036', ...
%     'F:\John B\Electrophys\Wild-Type\JBWT0024', ...
%     'F:\John B\Electrophys\Wild-Type\JBWT0025', ...
%     'F:\John B\Electrophys\Wild-Type\JBWT0029', ...
%     'F:\John B\Electrophys\Wild-Type\JBWT0032', ...
%     'F:\John B\Electrophys\Wild-Type\JBWT0033', ...
%     'D:\John B\Electrophysiology\Wild-Type Stim\JBWT0035' ...
    };
% sigs = {'visual', 'visual', 'rf', 'visual', 'visual', 'visual', 'visual', 'visual', 'visual'};
% recs = {5, [3 4], [2 4 7 9], [2 4], 4, [2 4], [4 5], [4 5], [3 4]};
% sigs = {'visual', 'visual', 'visual', 'visual'};
% recs = {[2 4], [4 5], [4 5], [3 4]};
% sigs = {'visual', 'visual'};
% recs = {4, 4};
sigs = repmat({'all'},1,5);
recs = {4 4 4 3 3};

assert(numel(dirs) == numel(sigs) && numel(sigs) == numel(recs));
currentDir = pwd;

load('analysis_20140116_seed','seed');

%%

for ii = 1:numel(dirs)
    rec = recs{ii};
    
    cd(dirs{ii});
    
    for jj = 1:numel(rec)
        try
            recordings = initRecordings;
%             analyseGollischMeisterResponses(recordings(rec(jj)),'forceclustered','yes','ignorenoise','yes','latencytype','peak','significance',sigs{ii});
            analyseJacobsNirenbergResponses(recordings(rec(jj)),'forceclustered','yes','ignorenoise','yes','latencytype','peak','significance',sigs{ii},'usesavedposteriors',false,'ifrfit','gauss','posteriorsonly',true,'seed',seed);
%             close all;
            
        catch err
            logMatlabError(err);
        end
    end
    
    cd(currentDir);
end

cd(currentDir);

%%

allRecordings = {};
for ii = 1:numel(dirs)
    rec = recs{ii};
    cd(dirs{ii});
    recordings = initRecordings;
    
    for jj = 1:numel(rec)
        recording = recordings(rec(jj));
        recording.parentDir = dirs{ii};
        allRecordings{end+1} = recording; %#ok<SAGROW>
    end
    
    cd(currentDir);
end

cd(currentDir);
    
analyseJacobsNirenbergResponses(allRecordings,'forceclustered','yes','ignorenoise','yes','latencytype','peak','significance',sigs{ii},'usesavedposteriors',true,'ifrfit','gauss','seed',seed);