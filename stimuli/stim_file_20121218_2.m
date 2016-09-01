repeats = 15;
speeds = 1:2;
angles = (0:11)*pi/6;
radius = 10;
centreX = 400;
centreY = 450;
maxT = 0;

nStimuli = repeats*length(angles);
getPixels = cell(2*nStimuli,1);
getExtraParams = cell(nStimuli,1);

stimulus = 0;
for ii = 1:repeats
    for kk = 1:12
        angle = angles(kk);

        stimulus = stimulus + 1;
        getExtraParams{stimulus} = struct( ...
            'angle',            angle,  ...
            'spatialPeriod',    16,     ...
            'temporalPeriod',   32);
        getPixels{2*stimulus-1} = @grating;
        getPixels{2*stimulus} = @grating;
    end
end

getExtraParams2 = getExtraParams(randperm(nStimuli));
getExtraParams3 = cell(2*nStimuli,1);

for ii = 1:nStimuli
    params = getExtraParams2{ii};
    maxT = maxT + 320;
    params.maxT = maxT;
    getExtraParams3{2*ii-1} = params;
    maxT = maxT + 60;
    params = struct('contrast',0,'luminance',0,'maxT',maxT);
    getExtraParams3{2*ii} = params;
end

stimulate(getPixels,[234 96 234+110 96+110],[],[],60,'C:\John Barrett\JBWT0009\Stim 3.mat',false,false,getExtraParams3);