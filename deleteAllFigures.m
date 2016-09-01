folders = {'Optogenetics' 'RD Pharma' 'Wild-Type Stim'};

currDir = pwd;
exptRegex = 'JB[AORW][DEGT][0-9]{4}';
stimRegex = '[0-9_]*(Stim|Spont) [0-9]+';

for ii = 1:3
    try
        cd(folders{ii});
        
        files = dir;
        fileNames = {files([files.isdir]).name};
        exptDirs = fileNames(~cellfun(@(s) isempty(regexp(s,exptRegex,'match')),fileNames));
        folder = pwd;
        
        for jj = 1:numel(exptDirs)
            try
                cd(exptDirs{jj});
                exptDir = pwd;
                
                files = dir;
                fileNames = {files([files.isdir]).name};
                stimDirs = fileNames(~cellfun(@(s) isempty(regexp(s,stimRegex,'match')),fileNames));
                
                for kk = 1:numel(stimDirs)
                    try
                        cd(stimDirs{kk});
                        
                        tic;
                        figs = dir('*.fig');
                        figs = {figs.name};
                        
                        if ~isempty(figs);
                            delete(figs{:});
                        end
                        
                        pngs = dir('*.png');
                        pngs = {pngs.name};
                        
                        if ~isempty(pngs)
                            delete(pngs{:});
                        end
                        
                        fprintf('Deleted %d figs and %d pngs in %f seconds\n',numel(figs),numel(pngs),toc);
                        
                        cd(exptDir);
                    catch err
                        logMatlabError(err);
                        cd(exptDir);
                    end
                end
                
                cd(folder);
            catch err
                logMatlabError(err);
                cd(folder);
            end
        end
        
        cd(currDir);
    catch err
        logMatlabError(err);
        cd(currDir);
    end
end