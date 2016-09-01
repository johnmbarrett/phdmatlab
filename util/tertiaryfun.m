function out = tertiaryfun(condition,thenfun,elsefun)
    if condition
        out = thenfun();
    else
        out = elsefun();
    end
end
        