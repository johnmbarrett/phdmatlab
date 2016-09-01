function getMoviePixels = movieStimulus(filename,nRepeats)
    if nargin < 1 || ~ischar(filename)
        [filename,filepath] = uigetfile({'*.mj2' 'Motion JPEG 2000'; '*.avi' 'AVI'; '.mpg' 'MPEG-1'; '.wmv' 'Windows Media Video'; '.asf' 'Windows Media Video'; '.asx' 'Windows Media Video'});
    else
        filepath = '';
    end
    
    if nargin < 2
        nRepeats = 1;
    end

    video = VideoReader([filepath filename]);
    nFrames = video.NumberOfFrames;

    getMoviePixels = @(T,varargin) getPixelsFromFrame(T,video,nFrames,nRepeats,struct());
end

function [pixels,extraParams] = getPixelsFromFrame(T,video,nFrames,nRepeats,extraParams)
    if T > nFrames*nRepeats
        pixels = NaN;
        return
    end
    
    tic;
    pixels = squeeze(read(video, mod(T-1,nFrames)+1));
    toc;
end
        