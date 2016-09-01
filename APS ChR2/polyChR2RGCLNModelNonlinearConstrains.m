function [c,ceq] = polyChR2RGCLNModelNonlinearConstrains(p)
    A = 3*p(1);
    B = 2*p(2);
    C = 1*p(3);
    
    roots = (-2*B+[1;-1]*sqrt(B^2-4*A*C))/(2*A);
    types = polyval([2*A B],roots);
    roots(types <= 0) = 0;
    
    c = [polyval(p(3:6),255)-1000;roots];
    c(~isreal(c)) = 0;
    ceq = [];
end