function createCStimBatchFile(periods,repeats,rgb,whiteTime,blackTime,dir,batchFile,outFilePrefix,exe)
    if nargin < 9
        exe = 'cstim.exe';
    end
    
    if nargin < 8
        outFilePrefix = 'test';
    end

    if nargin < 7
        batchFile = 'cstim.bat';
    end
    
    if nargin < 6
        dir = 'D:\backup\phd\cpp\cstim\cstim\x64\Debug';
    end
    
    if nargin < 5
        blackTime = 10;
    end
    
    if nargin < 4
        whiteTime = 5;
    end
    
    if nargin < 3
        rgb = [255 255 255];
    end
    
    if nargin < 2
        repeats = 10;
    end
    
    if nargin < 1
        periods = 2.^(4:6);
    end
    
    nTrials = numel(periods)*repeats;
    ps = repmat(periods',[repeats 1]);
    ps = ps(randperm(nTrials));
    
    fout = fopen(batchFile,'w');
    
    fprintf(fout,'%s\r\n',dir(1:2));
    fprintf(fout,'cd %s\r\n',dir(3:end));
    
    for ii = 1:nTrials
        fprintf(fout,'%s -r %d -g %d -b %d -p %d -t %d -o "%s_%d.txt"\r\n',exe,rgb(1),rgb(2),rgb(3),ps(ii),whiteTime,outFilePrefix,ii);
        fprintf(fout,'%s -r 0 -g 0 -b 0 -p 1000 -t %d -o dump.txt\r\n',exe,blackTime);
    end
    
    fprintf(fout,'del dump.txt\r\n');
    
    fclose(fout);
    
    disp(nTrials*(whiteTime+blackTime)/60);
end