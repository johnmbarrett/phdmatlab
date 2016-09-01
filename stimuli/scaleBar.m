function [pixels,extraParams,repeats] = scaleBar(T,x,y,window,extraParams) %#ok<INUSL>
    digits = cell(5,3);
    digits{1}  = [0 1 0; 1 0 1; 1 0 1; 1 0 1; 0 1 0];
    digits{2}  = [0 1 0; 1 1 0; 0 1 0; 0 1 0; 1 1 1];
    digits{3}  = [0 1 0; 1 0 1; 0 1 0; 1 0 0; 1 1 1];
    digits{4}  = [1 1 0; 0 0 1; 0 1 0; 0 0 1; 1 1 0];
    digits{5}  = [1 0 1; 1 0 1; 1 1 1; 0 0 1; 0 0 1];
    digits{6}  = [1 1 1; 1 0 0; 1 1 0; 0 0 1; 1 1 0];
    digits{7}  = [0 1 1; 1 0 0; 1 1 0; 1 0 1; 0 1 0];
    digits{8}  = [1 1 1; 0 0 1; 0 0 1; 0 1 0; 0 1 0];
    digits{9}  = [0 1 0; 1 0 1; 0 1 0; 1 0 1; 0 1 0];
    digits{10} = [0 1 0; 1 0 1; 0 1 1; 0 0 1; 1 1 0];
    
    if ~isfield(extraParams,'bgColour') || ~isnumeric(extraParams.colour) || ~ismember(numel(extraParams.colour),[1 3 4]) || any(extraParams.colour < 0 | extraParams.colour > 255)
        bgColour = 0;
    else
        bgColour = extraParams.bgColour;
    end
    
    if ~isfield(extraParams,'colour') || ~isnumeric(extraParams.colour) || ~ismember(numel(extraParams.colour),[1 3 4]) || any(extraParams.colour < 0 | extraParams.colour > 255)
        colour = 255;
    else
        colour = extraParams.colour;
    end
    
    if ~isfield(extraParams,'maxT') || ~isnumeric(extraParams.maxT) || numel(extraParams.maxT) ~= 1 || extraParams.maxT < 1 || isnan(extraParams.maxT)
        maxT = 1;
    else
        maxT = extraParams.maxT;
    end
    
    repeats = maxT;
    
    if numel(bgColour) ~= numel(colour)
        bgColour = 255-colour;
    end
    
    if T > maxT
        pixels = NaN;
        return
    end
    
    horizontal = isfield(extraParams,'horizontal') && all(extraParams.horizontal);
    if horizontal
        X = y;
        Y = x;
    else
        X = x;
        Y = y;
    end
    
    pix = zeros(Y,X);
    
    pix(:,end) = ones(Y,1);
    pix(mod(0:Y-1,2) == 0,end-(1:2)) = ones(ceil(Y/2),2);
    pix(mod(0:Y-1,10) == 0,end-(3:4)) = ones(ceil(Y/10),2);
    
    nGradations = ceil(Y/10);
    for ii = 0:nGradations-1
        rows = 10*(ii)+(1:5);
        rows = rows(rows < size(pix,1));
        
        zero = digits{1};
        
        pix(rows,end-(8:-1:6)) = zero(1:length(rows),:);
        
        nDigits = floor(log(ii)/log(10))+1;
        
        for jj = 1:nDigits
            digit = digits{mod(floor(ii/10^(jj-1)),10)+1};
            pix(rows,end-(8:-1:6)-4*jj) = digit(1:length(rows),:);
        end
    end
    
    if horizontal
        pix = pix';
        pix = flipud(pix);
    end
    
    pixels = squeeze(zeros([size(pix) numel(colour)]));
    for ii = 1:numel(colour)
        pixels(:,:,ii) = colour(ii)*pix + bgColour(ii)*(1-pix);
    end
end