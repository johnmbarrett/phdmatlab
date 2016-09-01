function [frameSamples,frameTimes] = extractPhotodiodeTimingsAPS(photodiodeFile,order,ripple,attenuation,lowPassFrequency,threshold,minFrameTime)
    % Extra stimulus frame times from the signal recorded from the
    % photodiode placed in front of the stimulation monitor
    %
    % PARAMATERS:
    %   photodiodeFile      - path to the file containing the photodiode data
    %   order               - order of the low-pass filter
    %   ripple              - passband ripple in dB
    %   attentuation        - stopband attentuation in dB
    %   lowPassFrequency    - cut-off frequency for low pass filter in Hz
    %   threshold           - threshold for determining stimuls timing
    %   minFrameTime        - minimum length of a stimulus frame in seconds
    %
    % OUTPUTS:
    %   frameSamples        - time of each frame in samples from the
    %                         beginning of the recording (1st sample = 1)
    %   frameTimes          - time of each frame in seconds from the
    %                         beginning of the recording
    
    % load data
    load(photodiodeFile,'Ch01_02','SamplingFrequency');
    x = Ch01_02;
    Fs = SamplingFrequency;
    
    % low pass filter
    [b,a] = ellip(order,ripple,attenuation,lowPassFrequency/(Fs/2));
    z = filtfilt(b,a,x);
    
    % find threshold-crossings
    frameSamples = find(diff(z > threshold));
    
    % remove crossings within minFrameTime of a previous crossing
    frameSamples = frameSamples([1 find(diff(frameSamples) > Fs*minFrameTime) + 1]);
    
    % convert to seconds
    frameTimes = frameSamples/Fs;
end