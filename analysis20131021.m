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
sigs = {'visual', 'visual', 'visual', 'rf', 'visual', 'visual', 'visual', 'visual', 'visual'};
recs = {[2 4 7 9], [2 4], 4, [2 4], [4 5], [4 5], [4 5], [3 4], [3 4]};
ages = [26 99 12 73 19 23 31 26 46];
lights = {[0 1 0 1],[0 1],1,[0 1],[0 0],[0 0],[0 0],[0 0],[0 0]};
sides = {'llrr','rr','u','rr','ll','ll','uu','ll','ll'};
tasks = {'ffff','ff','f','ff','fc','fc','fc','fc','fc'};
%%
cinfo = [];
linfo = [];
width = [];
contrast = [];
polarity = [];
age = [];
light = [];
side = [];
task = [];
retina = [];

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
            load(sprintf('%s\\%s_info',fileDir,recording.dataFile));
            
            n = size(unbiasedCountInfo,1);
            
            c = reshape(unbiasedCountInfo,n*4*2,1);
            l = reshape(unbiasedLatencyInfo,n*4*2,1);
            
            good = c <= 3 & l <= 3;
            
            if tasks{ii}(jj) == 'f'
                % this is wrong for one recording but eh
                w = reshape(repmat(2.^(3:6),[n 1 2]),size(c));
                k = ones(size(c));
            else
                w = 40*ones(size(c));
                k = reshape(repmat([1.0 0.4 0.2 0.1],[n 1 2]),size(c));
            end
            
            p = zeros(1,1,2);
            p(1,1,1) = 1;
            p(1,1,2) = 2;
            
            p = reshape(repmat(p,[n 4 1]),size(c));
            
            c = c(good);
            l = l(good);
            w = w(good);
            k = k(good);
            p = p(good);
            
            cinfo = [cinfo; c]; %#ok<AGROW>
            linfo = [linfo; l]; %#ok<AGROW>
            width = [width; w]; %#ok<AGROW>
            contrast = [contrast; k]; %#ok<AGROW>
            polarity = [polarity; p]; %#ok<AGROW>
            
            age = [age; ages(ii)*ones(size(c))]; %#ok<AGROW>
            light = [light; lights{ii}(jj)*ones(size(c))]; %#ok<AGROW>
            side = [side; repmat(sides{ii}(jj),size(c))]; %#ok<AGROW>
            task = [task; repmat(tasks{ii}(jj),size(c))]; %#ok<AGROW>
            retina = [retina; r*ones(size(c))]; %#ok<AGROW>
        end
        
        r = r + 1;
        cd ..;
    catch err
        logMatlabError(err);
        cd(currentDir);
    end
end

dinfo = cinfo-linfo;

%%
save('aggregate_info','dinfo','cinfo','linfo','width','contrast','polarity','age','light','side','task','retina');