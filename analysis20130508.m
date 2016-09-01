dirs = {...'../JBWT0028', '../JBWT0029', '../JBWT0026',
    'F:\John B\Electrophys\JBWT0024', 'F:\John B\Electrophys\JBWT0025'};
    %'C:\temp\retina data\JBWT0027'};
sigs = {...'visual' 'visual', 'visual', 
    'rf', 'visual'};
%     'visual'};
recs = {...[2 4], [2 4], 4, 
    [2 4 7 9], [2 4]};
    %[7 9]};

assert(numel(dirs) == numel(sigs) && numel(sigs) == numel(recs));
currentDir = pwd;

for ii = 1:numel(dirs)
    rec = recs{ii};
    
    cd(dirs{ii});
    
    for jj = 1:numel(rec)
        try
            recordings = initRecordings;
            analyseJacobsNirenbergResponses(recordings(rec(jj)),'forceclustered','yes','ignorenoise','yes','significance',sigs{ii});
            close all;
        catch err
            logMatlabError(err);
            cd(currentDir);
        end
    end
end

cd(currentDir);