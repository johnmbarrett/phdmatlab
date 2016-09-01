function [histograms,edgess,repeats] = rastersToPSRH(rasters,stimulusLengths,bw)
    histograms = cell(size(rasters));
    edgess = cell(size(rasters));
    repeats = zeros(size(rasters));
    
    for ii = 1:numel(rasters)
        raster = rasters{ii};
        nLines = numel(raster);
        
        subscript = cell(1,ndims(repeats));
        [subscript{:}] = ind2sub(size(repeats),ii);
        repeats(ii) = nLines;
        
        stimulusLengthSubscript = subscript(2:end);
        stimulusLength = stimulusLengths(sub2ind(size(stimulusLengths),stimulusLengthSubscript{:}));
        stimulusLengthMS = ceil((stimulusLength+1)*10);
        edges = linspace(-1,stimulusLengthMS*bw,stimulusLengthMS+1); % 100ms bins
        edgess{ii} = edges;
        
        for jj = 1:nLines
            histogram = histc(raster{jj},edges);
            
            if isempty(histogram)
                continue;
            end
            
            histogram = reshape(histogram,numel(histogram),1)/nLines;
            
            if isempty(histograms{ii})
                histograms{ii} = histogram;
            else
                histograms{ii} = histograms{ii} + histogram;
            end
        end
    end
end