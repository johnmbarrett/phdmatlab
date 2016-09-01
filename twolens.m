function [L, dj, di, table] = twolens(do, f1, f2, M)
    di = (f1*do)./(do-f1);
    L  = (di*f2)./(do*M) + di + f2;
    dj = ((L-di)*f2)./(L-di-f2);
   
    if isrow(do)
        table = [do; L; dj; di; do+L+dj];
    elseif iscolumn(do)
        table = [do L dj di do+L+dj];
    end
end

% function twolens()
%     % constants
%     M = 0.2;
%     f1 = 100;
%     f2 = 50;
% 
%     % variables
%     do = 0;
%     dj = 0;
%     k = 0.5;
%     L = 0;
% end
% 
% % equations
% 
% function [k, dj] = moveObjectOr1stLens(do, L, M, f1)
%     k = f1*do/L(do - f1);
%     dj = (1-k)*M*do/k;
% end
% 
% function [do, dj, L, k] = move1stLens(delta, do, L, M, f1)
%     do = do+delta;
%     L = L-delta;
%     [k, dj] = moveObjectOr1stLens(do, L, M, f1);
% end
% 
% function [k, dj] = moveObject(delta, do, L, M, f1)
%     [k, dj] = moveObjectOr1stLens(do+delta, L, M, f1);
% end
% 
% function [dj, L, k] = move2ndLens(delta, L, k, f2)
%     k = k*L/(L+delta);
%     L = L+delta;
%     dj = (1-k)*L*f2/((1-k)*L-f2);
% end
% 
% function [do, dj] = adjustK(k, L, f1, f2, M)
%     do = k*L*f1/(k*L-f1);
