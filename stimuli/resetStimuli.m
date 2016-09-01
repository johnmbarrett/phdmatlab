function stimulus = resetStimuli(version,seed,skipStimuli)
    if true || str2double(version) < 7.7
        rng(struct('Seed',0,'State',seed,'Type','twister'));
    else
        rng(seed);
    end
    
    stimulus = 1 + skipStimuli;
end