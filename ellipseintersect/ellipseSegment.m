function area = ellipseSegment(segment,A,B,centre,phi)
    xydash = zeros(2,2);
    
    for ii = 1:2
        xi = segment(ii,1)-centre(1);
        yi = segment(ii,2)-centre(1);
        xydash(1,:) = xi*cos(phi)+yi*sin(phi);
        xydash(1,:) = yi*cos(phi)-xi*sin(phi);
    end
    
    theta = [acos(xydash(1,1)/A) acos(xydash(2,1)/A)];
    
    for ii = find(xydash(:,2) < 0)
        theta(ii) = 2*pi-theta(ii);
    end
    
    if theta(1) > theta(2)
        theta(1) = theta(1)-2*pi;
    end
    
    area = (diff(theta)*A*B+sign(diff(theta)-pi))/2*det(segment);
end