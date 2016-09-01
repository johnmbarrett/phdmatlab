function Z = gauss3d(X,Y,T,sdx,sdy,sdt,phi,normalise)
    a = cos(phi)^2/(2*sdx^2) + sin(phi)^2/(2*sdy^2);
    b = sin(2*phi)/(4*sdy^2) - sin(2*phi)/(4*sdx);
    c = sin(phi)^2/(2*sdx^2) + cos(phi)^2/(2*sdy^2);
    Zs = exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
    Zt = exp(-T.^2/(2*sdt^2));
    Z = repmat(Zs,[1 1 numel(Zt)]).*repmat(reshape(Zt,1,1,numel(Zt)),size(Zs,1),size(Zs,2));
    
    if nargin > 7 && normalise
        Z = Z/max(max(max(Z)));
    end
end