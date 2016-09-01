function Z = plogp(X,Y)
    if nargin < 2
        Y = X;
    end
    
    % TODO : check all past usages of this function to see if they were
    % assuming nats or bits
    Z = X.*log(Y);
    Z(isnan(Z)) = 0;
end