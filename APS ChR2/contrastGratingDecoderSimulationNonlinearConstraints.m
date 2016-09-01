function [c,ceq] = contrastGratingDecoderSimulationNonlinearConstraints(p)
    c = 500*gampdf(80,p(3),p(4))-1;
    ceq = [];
end