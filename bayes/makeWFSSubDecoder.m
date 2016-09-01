function subDecoder = makeWFSSubDecoder(decoder,featureIndices,varargin)
    subDecoder = struct('latencyMatrix',decoder.latencyMatrix(featureIndices,featureIndices,:));
end