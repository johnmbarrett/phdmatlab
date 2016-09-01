function P = generateLegendrePolynomials(upto,outfile)
    assert(upto >= 0,'First input argument must be non-negative.');
    
    P = cell(upto+1,1);
    P{1} = 1;
    
    if upto > 0
        P{2} = [1 0];
    end
    
    for ii = 2:upto
        P{ii+1} = ([(2*ii-1)*P{ii} 0]-[0 0 (ii-1)*P{ii-1}])/ii;
    end
end
        
    
    