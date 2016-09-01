function Z = gauss2d(X,Y,sdx,sdy,phi,normalise)
    a = cos(phi)^2/(2*sdx^2) + sin(phi)^2/(2*sdy^2);
    b = sin(2*phi)/(4*sdy^2) - sin(2*phi)/(4*sdx^2);
    c = sin(phi)^2/(2*sdx^2) + cos(phi)^2/(2*sdy^2);
    Z = exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
    
    if nargin > 5 
        switch normalise
            case 0
                return;
            case 1
                Z = Z/max(max(Z));
            case 2
                Z = Z/sum(Z(:));
        end
    end
end