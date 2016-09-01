function pr = paethPredictor(a,b,c)
    assert(isequal(size(a),size(b) && isequal(size(b),size(c))));

    p = a+b-c;
    pa = abs(p-a);
    pb = abs(p-b);
    pc = abs(p-c);
    pr = nan(size(p));
    
    pr(pa <= pb & pa <= pc) = a;
    pr(isnan(pr) & pb <= pc) = b;
    pr(isnan(pc)) = c;
end
    