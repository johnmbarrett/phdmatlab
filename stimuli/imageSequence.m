function f = imageSequence(imagePaths,sequence)
    if nargin ~= 2
        error 'You must provide a list of images and a presentation order';
    end
    
    images = cell(size(imagePaths));
    
    for ii = 1:numel(imagePaths)
        image = squeeze(double(imread(imagePaths{ii})));
        
        if isequal(unique(image),[0;1])
            image = 255*image;
        end
        
        images{ii} = image;
    end
    
    sequence(:,2) = cumsum(sequence(:,2));
    
    extraParams = struct([]);
    extraParams(1).sequence = sequence;
    extraParams(1).images = images;
    extraParams(1).maxT = sequence(end,2);
    
    f = @(T,X,Y,W,P) getImageSequence(T,X,Y,W,extraParams);
end
    
function [pixels,extraParams,repeats,markerIndex] = getImageSequence(T,~,~,~,extraParams)
    if T > extraParams.maxT
        markerIndex = NaN;
        repeats = 0;
        pixels = NaN;
        return
    end

    nextImage = find(extraParams.sequence(:,2) >= T,1);
    
    if nextImage == 1
        repeats = extraParams.sequence(1,2);
    else
        repeats = diff(extraParams.sequence(nextImage-[1 0],2));
    end
    
    pixels = extraParams.images{extraParams.sequence(nextImage,1)};
    
    markerIndex = mod(nextImage-1,2)+1;
end