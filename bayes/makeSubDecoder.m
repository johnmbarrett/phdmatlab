function newValues = makeSubDecoder(oldValues,featureIndices,nFeatures)
    s = size(oldValues);
    
    d = find(s == nFeatures,1);
    
    if isempty(d)
        newValues = oldValues;
        return;
    end
    
    coords = repmat({':'},1,ndims(oldValues));
    coords{d} = featureIndices;
    
    newValues = oldValues(coords{:});
end