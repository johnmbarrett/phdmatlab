function [c,ceq] = contrastGratingsMutualInformationNonlinearConstraint(p)
    c = -polyval(p(3:6),1:6);
    ceq = [];
end