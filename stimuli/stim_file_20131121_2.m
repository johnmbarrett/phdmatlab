X = 334-(140-56)/2;
Y = 247-(140-56)/2;
repeats = 150;
fps = 30;
stimulusDir = 'D:\John B\Stimuli\';
markerRect = [0 0 120 640];
saveFile = 'D:\John B\Electrophysiology\Wild-Type Stim\JBWT0040\Stim 4.mat';

textureRect = [X Y X+140 Y+140];

time = 0;

grey = [stimulusDir 'grey.png'];
singleImages = dir([stimulusDir 'geoffrey\all_12_images\*.png']);
images = cell(13,1);

for ii = 1:numel(singleImages)
    images{ii} = [stimulusDir 'geoffrey\all_12_images\' singleImages(ii).name];
end

images{13} = grey;

sequence = zeros(12*2*repeats,2);

for ii = 1:repeats
    sequence((ii-1)*24+1:2:ii*24,1) = randperm(12);
    sequence((ii-1)*24+1:2:ii*24,2) = 1;
    sequence((ii-1)*24+2:2:ii*24,1) = 13;
    sequence((ii-1)*24+2:2:ii*24,2) = 23;
end

time = time + sum(sequence(:,2));

disp(time/(fps*60));

shortExposure = imageSequence(images,sequence);

sequence = zeros(12*2*repeats,2);

for ii = 1:repeats
    sequence((ii-1)*24+1:2:ii*24,1) = randperm(12);
    sequence((ii-1)*24+1:2:ii*24,2) = 12;
    sequence((ii-1)*24+2:2:ii*24,1) = 12;
    sequence((ii-1)*24+2:2:ii*24,2) = 12;
end

time = time + sum(sequence(:,2));

disp(time/(fps*60));

longExposure = imageSequence(images,sequence);

subfolders = { ...
    'microsaccade-grating-100' ...
    'microsaccade-grating-200' ...
    'microsaccade-grating-800' ...
    'microsaccade-imageNat-01' ...
    'microsaccade-imageNat-02' ...
    'microsaccade-imageNat-03' ...
    'microsaccade-imageNat-sketch-01' ...
    'microsaccade-imageNat-sketch-02' ...
    'microsaccade-imageNat-sketch-03' ...
    'microsaccade-texture-01' ...
    'microsaccade-texture-02' ...
    'microsaccade-texture-03' ...
    };

images = cell(12*12+1,1);

for ii = 1:12
    subfolder = subfolders{ii};
    
    for jj = 1:12
        images{(ii-1)*12+jj} = [stimulusDir 'geoffrey\' subfolder sprintf('\\shifted_%02d.png',jj-1)];
    end
end

images{145} = grey;

sequence = zeros(13*12*repeats,2);

for ii = 1:repeats
    stimulusOrder = 12*(randperm(12)-1)';
    stimulusOrder = kron(stimulusOrder,ones(13,1))+repmat((1:13)',12,1);
    stimulusOrder(13:13:end) = 145;
    sequence(13*12*(ii-1)+1:13*12*ii,1) = stimulusOrder;
end

sequence(:,2) = repmat([ones(12,1); 12],12*repeats,1);

time = time + sum(sequence(:,2));

disp(time/(fps*60));

shortSaccades = imageSequence(images,sequence);

sequence = [zeros(61*12,1) ones(61*12,1)];

sequence(61:61:end,1) = 145;
sequence(61:61:end,2) = 60;

stimulusOrder = randperm(12);

for ii = 1:12
    sequence(61*(ii-1)+1:61*ii-1,1) = 12*(stimulusOrder(ii)-1)+repmat((1:12)',5,1);
end

time = time + sum(sequence(:,2));

longSaccades = imageSequence(images,sequence);

disp(time/(fps*60));

stimulate({shortExposure; longExposure; shortSaccades; longSaccades},textureRect,[],markerRect,255,0,fps,saveFile,false,false,false,cell(4,1));