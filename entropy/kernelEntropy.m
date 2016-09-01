function h = kernelEntropy(x)
   n = numel(x);
   h = 0;

   options = optimset(optimset(@fmincon),'Algorithm','active-set');
   
   for jj = 1:n
       Ihat = @(w) -sum(plogp(sum(tpdf((x(setdiff(1:n,jj))-x(jj))/w,w))/((n-1)*w)))/n;

       w = fmincon(Ihat,sqrt(n),[],[],[],[],1,n,[],options);

       Ihatw = Ihat(w);
       h = h + Ihatw/n;
   end
end