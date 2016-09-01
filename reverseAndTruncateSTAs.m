recordings = initRecordings;

for recording = recordings
    cd(getAnalysisOutputDir(recording));
    gifs = dir('*.gif');
    gifs = {gifs.name};
    
    for jj = 1:length(gifs)
        gif = gifs{jj};
        I = imread(gif,'gif','frames','all');
        
        if size(I,4) == 90
            if ~movefile(gif,regexprep(gif,'.gif','.old'),'f')
                goping;
                error('Bums');
            end
            
            I = I(:,:,:,30:-1:1);
            imwrite(I,gif,'LoopCount',1,'DelayTime',0.1);
        end
    end
    
    cd ..;
end

goping;