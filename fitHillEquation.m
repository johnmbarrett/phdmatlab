function [threshold,exponent,maxY,R,J,CovB,MSE] = fitHillEquation(X,Y)
    assert(size(X,1) == size(Y,1),'X and Y must have the same number of rows');
    assert(ismember(size(X,2),[1 size(Y,2)]),'X must have one value for either every row of or every value in Y');
    
    maxY = median(Y(X == max(max(X)),:));
    
    if size(X,2) == 1
        X = repmat(X(:),size(Y,2),1);
    else
        X = reshape(X,numel(X),1);
    end
    
    Y = reshape(Y,numel(Y),1);
    
    [p,R,J,CovB,MSE] = nlinfit(X,Y,@(p,x) maxY.*x.^p(1)./(p(2).^p(1) + x.^p(1)),[1 maxY/2]);
    exponent = p(1);
    threshold = p(2);
end