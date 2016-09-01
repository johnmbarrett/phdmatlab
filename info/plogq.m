function z = plogq(x,y)
    z = x.*log(y);
    z(isnan(z)) = 0;
end