screenRect = [1920 0 1920+800 600];

%%

saveFile1 = 'D:\John B\Electrophysiology\Optogenetics\JBOG0045\Stim 9.mat';
getPixels = colourSequence({255;0},[0 0 800 600],120,7200,false,0);
stimulate(getPixels,[0 0 800 600],screenRect,[],0,0,60,saveFile1,false,false,false,struct([]));

%%

durations = 2.^(0:6);
intensities = [0 31:32:255];

on = cell(7,9);
off = cell(7,1);

for ii = 1:7
    t = durations(ii);
    T = 120-durations(ii);
    off{ii} = colourSequence(0,[0 0 800 600],T,T,false,0);
    
    for jj = 1:9
        on{ii,jj} = colourSequence(intensities(jj),[0 0 800 600],t,t,false,0);
    end
end

%%

reps = 10;

getPixels = cell(7*9*reps*2,1);

for ii = 1:reps
    indices = randperm(7*9);
    
    for jj = 1:(7*9)
        [y,x] = ind2sub([7 9],indices(jj));
        getPixels(10*(ii-1)+2*jj-1) = on(y,x);
        getPixels(10*(ii-1)+2*jj) = off(y);
    end
end

%%

saveFile2 = NaN;

stimulate(getPixels,[0 0 800 600],screenRect,[],0,0,60,saveFile2,false,false,false,cell(size(getPixels)));