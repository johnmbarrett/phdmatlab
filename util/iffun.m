function out = iffun(condition,thenfun,elsefun)
    if condition
        out = thenfun();
    else
        out = elsefun();
    end
end
        