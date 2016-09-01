function [Z,X,Y] = dog(X,Y,sdXon,sdYon,sdXoff,sdYoff,phi,k,normalise,offset)
    if nargin < 10
        offset = 0;
    end
    
    if nargin < 9
        normalise = false;
    end
    
    if nargin < 8
        k = 0.5;
    end
    
    if nargin < 7
        phi = 0;
    end
    
    if nargin < 6
        sdYoff = 2;
    end
    
    if nargin < 5
        sdXoff = 2;
    end
    
    if nargin < 4
        sdYon = 1;
    elseif isempty(sdYoff) || isnan(sdYoff)
        % can specify surround standard deviation separately or as a fixed
        % multiple of centre standard deviation
        sdYoff = sdXoff*sdYon;
        sdXoff = sdXoff*sdXon;
    end
    
    if nargin < 3
        sdXon = 1;
    end
    
    if sdXon >= sdXoff
        temp = sdXon;
        sdXon = sdXoff;
        sdXoff = temp;
    end
    
    if sdYon >= sdYoff
        temp = sdYon;
        sdYon = sdYoff;
        sdYoff = temp;
    end
    
    if nargin == 1
        error('You must provide both X and Y meshers or neither');
    end
    
    if nargin < 1
        [Y,X] = ndgrid(-5*sdXoff:5*sdXoff,-5*sdYoff:5*sdYoff);
    end
    
    Z = gauss2d(X,Y,sdXon,sdYon,phi,true)-k*gauss2d(X,Y,sdXoff,sdYoff,phi,true)+offset;
    
    if normalise && max(max(Z)) ~= 0
        Z = Z/max(max(Z));
    end
end