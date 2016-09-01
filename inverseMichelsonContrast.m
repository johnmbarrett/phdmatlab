function [Imax,Imin] = inverseMichelsonContrast(contrast,luminance)
    if nargin < 2
        luminance = 127.5;
    end
    
    if nargin < 1
        contrast = 1;
    end
    
    Imax = (1 + contrast)*luminance;
    Imin = (1 - contrast)*luminance;
end