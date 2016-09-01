function position = ellipseCase(A,B,AA,BB,CC,DD,EE,FF)
    d = AA(1)*(CC(1)*FF(1)-EE(1)^2) - ...
        (CC(1)*DD(1)^2-2*BB(1)*DD(1)*EE(1)+FF(1)*BB(1)^2);
    
    b = (AA(1)*(CC(2)*FF(2)-EE(2)^2) + ...
        2*BB(1)*(EE(2)*DD(2)-FF(2)*BB(2)) + ...
        2*DD(1)*(EE(2)*BB(2)-CC(2)*DD(2)) + ...
        CC(1)*(AA(2)*FF(2)-DD(2)^2) + ...
        2*EE(1)*(BB(2)*BB(2)-AA(2)*EE(2)) + ...
        FF(1)*(AA(2)*CC(2)-BB(2)^2))/d;
    
    a = (AA(1)*(CC(1)*FF(2)-2*prod(EE)+FF(1)*CC(2)) + ...
        2*BB(1)*(EE(1)*DD(2)-FF(1)*BB(2)+DD(1)*EE(2)) + ...
        2*DD(1)*(EE(1)*BB(2)-CC(1)*DD(2)) - ...
        (BB(1)^2*FF(2)+DD(1)^2*CC(2)+EE(1)^2*B(1,1)) + ...
        CC(1)*FF(1)*AA(2))/d;
    
    c = (AA(2)*(CC(2)*FF(2)-EE(2)^2)-(BB(2)^2*FF(2)-2*BB(2)*DD(2)*EE(2)+DD(2)^2*CC(2)))/d;
    
    s4 = -27*c^3+18*c*a*b+a^2*b^2-4*a^3*c-4*b^3;
    
    if s4 < 0
        position = 2;
        return;
    end
    
    s1 = a;
    s2 = a^2-3*b;
    s3 = 3*a*c+b*a^2-4*b^2;
    
    if s4 > 0
        if ~(s1 > 0 && s2 > 0 && s3 > 0)
            position = 3;
            return;
        end
        
        u = (-a-sqrt(s2))/3;
        v = (-a+sqrt(s2))/3;
        
        M = u*A+B;
        N = v*A+B;
        
        if (M(2,2)*det(M) > 0 && det(M(2:3,2:3)) > 0) ...
        || (N(2,2)*det(N) > 0 && det(N(2:3,2:3)) > 0)
            position = 4;
        else
            position = 1;
        end
        
        return;
    end
    
    if s1 > 0 && s2 > 0 && s3 < 0
        position = 6;
        return;
    end
    
    if s1 > 0 && s2 > 0 && s3 > 0
        beta = (9*c-a*b)/(2*s2);
        alpha = (4*a*b-a^3-9*c)/s2;
        
        M = beta*A+B;
        N = alpha*A+B;
        
        if det(M(1:2,1:2)) <= 0
            position = 8;
            return;
        end
        
        if det(M(2:3,2:3))+det(M([1 3],[1 3])) > 0
            position = 4;
            return;
        end
        
        if det(N(1:2,1:2)) > 0
            position = 5;
        else
            position = 7;
        end
        
        return;
    end
    
    alpha = -a/3;
    M = alpha*A+B;
    
    % TODO : in the original paper there's a condition of the form P|(P&Q),
    % which reduces to P, which is almost certainly a typo.  Check with the
    % paper it references.
    if det(M) ~= 0 && det(M(1:2,1:2)) <= 0
        position = 10;
        return;
    end
    
    if det(M) ~= 0 && det(M(1:2,1:2)) > 0
        position = 9;
        return;
    end
    
    position = 7;
end