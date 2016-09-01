function y = plogp(x)
    y = x.*log(x);
    y(isnan(y)) = 0;
end