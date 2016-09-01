function Z = gauss1d(X,sd,normalise)
    Z = exp(-X/sd);
    
    if nargin > 2 
        switch normalise
            case 0
                return;
            case 1
                Z = Z/max(Z);
            case 2
                Z = Z/sum(Z);
        end
    end
end