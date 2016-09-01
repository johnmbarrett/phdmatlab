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
            'maxT',             320,    ...
            'spatialPeriod',    16,     ...
            'temporalPeriod',   32);
        getPixels{2*stimulus-1} = @grating;
        getPixels{2*stimulus} = @grating;
    end
end

getExtraParams2 = getExtraParams(randperm(nStimuli));
getExtraParams3 = cell(2*nStimuli,1);

for ii = 1:nStimuli
    getExtraParams3{2*ii-1} = getExtraParams2{ii};
    params = struct('contrast',0,'luminance',0,'maxT',60);
    getExtraParams3{2*ii} = params;
end

stimulate(getPixels,[194 67 194+222 67+222],[],[],60,'C:\John Barrett\JBWT0010\Stim 3.mat',false,false,getExtraParams3);