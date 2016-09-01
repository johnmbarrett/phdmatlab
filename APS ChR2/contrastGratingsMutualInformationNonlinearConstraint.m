function [c,ceq] = contrastGratingsMutualInformationNonlinearConstraint(p)
    c = (-2*p(4)+[-1 1]*sqrt(4*p(4)-12*p(3)*p(5)))/(2*p(3));
    c(~isreal(c)) = 0;
    ceq = [];
end