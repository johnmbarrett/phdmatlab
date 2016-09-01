function [Vt,t] = gjsim(V0,Cm,g,t2,t1)
    c = 1e-2;
    m = 1e-3;
    u = 1e-6;
    
    if nargin < 5
        t1 = 0;
    end
    
    if nargin < 4
        t2 = 1;
    end

    if nargin < 3
        g = 0.05*m/c^2;
    end

    if nargin < 2
        Cm = 1*u/c^2;
    end

    k = g/Cm;

    if nargin < 1
        V0 = [-60;0;-20]*m;
    end;
    
    function dV = gapJunction(~,V)
        dV = zeros(size(V));
        
        for ii = 1:numel(V)
            for jj = 1:numel(V);
                dV(ii) = dV(ii) + k*(V(jj)-V(ii));
            end
        end
    end
    
    [t,Vt] = ode45(@gapJunction,[t1 t2],V0);
end

