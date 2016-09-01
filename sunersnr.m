function [SNR, A, SDnoise, Wbar] = sunersnr(W)
%SUNERSNR Calculates SNR using Suner's forumla
%   Calculates SNR using the formula defined in Suner et al (2005), IEEE
%   Trans Neur Sys Rehab Eng.  Input W is a k*N matrix of k spikes
%   waveforms, each comprising N samples.  Outputs are the signal-to-noise
%   ratio SNR, peak-to-peak voltage A, standard deviation SDnoise and mean
%   spike waveform Wbar

    Wbar = mean(W);
    A = max(Wbar)-min(Wbar);
    err = W - ones(size(W,1),1)*Wbar;
    SDnoise = std(reshape(err,numel(err),1));
    SNR = A/(2*SDnoise);

end

