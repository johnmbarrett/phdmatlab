function [a,AA,BB,CC,DD,EE,FF] = ellipseImplicitPolynomial(A,B,x0,y0,phi)
    n = numel(A);
    
    a = zeros(3,3,n);
    
    A = A.^2;
    B = B.^2;
    
    sp = sin(phi);
    cp = cos(phi);
    
    a(1,1,:) = cp.^2./A + sp.^2./B;
    a(1,2,:) = sp.*cp./A - sp.*cp./B;
    a(2,2,:) = sp.^2./A + cp.^2./B;
    
    a(1,3,:) = sp.*(x0.*cp-y0.*sp)./B - cp.*(x0.*cp+y0.*sp)./A;
    a(2,3,:) = cp.*(x0.*sp-y0.*cp)./B - sp.*(x0.*cp+y0.*sp)./A;
    a(3,3,:) = (x0.*cp+y0.*sp).^2./A + (x0.*sp-y0.*cp).^2./B - 1;
    
    AA = squeeze(a(1,1,:));
    BB = squeeze(2*a(1,2,:));
    CC = squeeze(a(2,2,:));
    DD = squeeze(2*a(1,3,:));
    EE = squeeze(2*a(2,3,:));
    FF = squeeze(a(3,3,:));
    
    a(2,1,:) = a(1,2,:);
    a(3,1,:) = a(1,3,:);
    a(3,2,:) = a(2,3,:);
end