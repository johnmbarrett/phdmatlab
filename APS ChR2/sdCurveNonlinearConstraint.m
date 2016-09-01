function [c,ceq] = sdCurveNonlinearConstraint(params,threshold)
    ceq = [];
    
    c = zeros(2,1);
    
    c(1) = -params(4)/params(3);
    
    s = threshold/params(1)+0.5;
    k = log(s/(1-s));
    
    c(2) = -k/(params(2)*params(3));
end 