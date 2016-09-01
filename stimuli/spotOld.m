function f = spot(centre,radius,colour,bgColour,maxT)
    if nargin < 5
        maxT = Inf;
    end
    
    if nargin < 4
        bgColour = 0;
    end
    
    if nargin < 3
        colour = 255;
    end
    
    colourDepth = numel(colour);
    
    if colourDepth ~= numel(bgColour)
        bgColour = zeros(size(colour));
    end
    
    if nargin < 2
        radius = 0.5;
    end
    
    if nargin < 1
        error 'You must specify a centre location';
    end
    
    if numel(centre) ~= 2
        error 'You must specify the centre location as an (x,y) co-ordinate pair';
    end
    
%     circle = Circle(radius);
    space = linspace(-radius,radius,ceil(radius*2));
    [x,y] = meshgrid(space,space);
    circle = (x.^2 + y.^2) <= radius;
    
    circlePixels = zeros([size(circle) colourDepth]);
    
    for ii = 1:colourDepth
        circlePixels(:,:,ii) = circle*colour(ii) + (1-circle)*bgColour(ii);
    end
    
    extraParams = struct('colourDepth',colourDepth,'maxT',maxT);
    extraParams.bgColour = bgColour;
    extraParams.centre = centre;
    extraParams.circlePixels = circlePixels;
    
    f = @(T,X,Y,W,P) getSpot(T,X,Y,W,extraParams);
end

function [pixels,extraParams] = getSpot(T,X,Y,window,extraParams)
    if T > extraParams.maxT
        pixels = NaN;
        return;
    end
    
    pixels = ones(X,Y,extraParams.colourDepth);
    
    centre = extraParams.centre;
	circle = extraParams.circlePixels;
    circleSize = size(circle);
    pixelSize = size(pixels);
    
    cxMin = 1;
    cxMax = circleSize(1);
    cyMin = 1;
    cyMax = circleSize(2);
    
    spill = false;
    
    pxMin = centre(1)-ceil(circleSize(1)/2)+1;
    
    if pxMin < 1
        spill = true;
        
        cxMin = 2-pxMin;
        pxMin = 1;
    end
    
    pxMax = centre(1)+floor(circleSize(1)/2);
    
    if pxMax > pixelSize(1)
        spill = true;
        
        cxMax = cxMax - pxMax + pixelSize(1);
        pxMax = pixelSize(1);
    end
    
    pyMin = centre(2)-ceil(circleSize(2)/2)+1;
    
    if pyMin < 1
        spill = true;
        
        cyMin = 2-pyMin;
        pyMin = 1;
    end
    
    pyMax = centre(2)+floor(circleSize(2)/2);
    
    if pyMax > pixelSize(1)
        spill = true;
        
        cyMax = cyMax - pyMax + pixelSize(2);
        pyMax = pixelSize(2);
    end
    
    if spill
        fprintf('Warning: spot (size [%d %d] centre [%d %d]) extends beyonds the borders of the texture rectangle (size [%d %d]), parts of the spot will be cropped to fit\n',circleSize(1),circleSize(2),centre(1),centre(2),pixelSize(1),pixelSize(2));
    end
    
    for ii = 1:extraParams.colourDepth
        pixels(:,:,ii) = extraParams.bgColour(ii)*pixels(:,:,ii);
        pixels(pxMin:pxMax,pyMin:pyMax,ii) = circle(cxMin:cxMax,cyMin:cyMax,ii);
    end
end
    