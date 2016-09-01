function gifToPng(newDir)
    if nargin < 1
        newDir = pwd;
    end
    
    oldDir = pwd;
    
    clean = onCleanup(@() cd(oldDir));
    
    cd(newDir);
    
    gifs = dir('*.gif');
    names = {gifs.name};

    for ii = 1:numel(names)
        [~,name] = fileparts(names{ii});
        I = imread([name '.gif']);
        imwrite(I,[name '.png'],'png');
        delete(names{ii});
    end
end