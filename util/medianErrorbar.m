function [h,m,l,u] = medianErrorbar(X,Y,dim,varargin)
    if nargin < 2
        Y = X;
        X = 1:size(Y,2);
    end
    
    sy = size(Y);
    
    if nargin < 3
        dim = [];
    elseif ischar(dim)
        varargin = [{dim} varargin];
        dim = [];
    end
    
    sx = size(X);
    
    if isvector(X)
        xdim = setdiff(find(sy == numel(X)),dim);
        
        if isempty(xdim)
            error('X must be a vector of length equal to one of the dimensions of Y or else there must be one X data point for each condition.')
        end
        
        if isempty(dim)
            dims = setdiff(1:ndims(Y),xdim);
            dim = dims(1);
        end
        
        Y = permute(Y,[dim xdim setdiff(1:ndims(Y),[dim xdim])]);
        sy = size(Y);
        X = repmat(X(:),[1 sy(3:end)]);
    else
        assert(isequal(sx,sy(2:end)),'X must be a vector of length equal to one of the dimensions of Y or else there must be one X data point for each condition.');
    end
    
    if numel(sy) == 2
        sy(3) = 1;
    end
    
    q = prctile(Y,[50 25 75]);
    q = reshape(q,[3 prod(sy(2:end))]);
    m = reshape(squeeze(q(1,:)),sy(2:end));
    l = reshape(squeeze(q(1,:)-q(2,:)),sy(2:end));
    u = reshape(squeeze(q(3,:)-q(1,:)),sy(2:end));
    
    h = errorbar(X,m,l,u,varargin{:});
end