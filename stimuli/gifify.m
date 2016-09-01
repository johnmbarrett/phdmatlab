function Y = gifify(X)
    if ndims(X) > 3
        error '3D animated gifs are not supported at this time.';
    end
    
    Y = reshape(X,size(X,1),size(X,2),1,size(X,3));
end