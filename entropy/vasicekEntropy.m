function Vmn = vasicekEntropy(x,m,unbiased,base)
   if nargin < 4
       base = exp(1);
   end

   n = numel(x);
   y = sort(x);
   
   y1 = [y(1)*ones(m,1); y];
   y2 = [y; y(n)*ones(m,1)];
   
   z = y2((m+1):end)-y1(1:n);
   
   Vmn = sum(log(n*z/(2*m)))/n;
   
   if nargin > 2 && unbiased
       Vmn = Vmn - log(n) + log(2*m) - (1-2*m/n)*psi(2*m) + psi(n+1) - 2*sum(psi((1:m)+m-1))/n;
   end  
   
   Vmn = Vmn/log(base);
end