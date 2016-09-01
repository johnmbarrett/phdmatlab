function y = safelog(x)
    y = x;
    y(x <= 0) = realmin;
    y = log(y);
end