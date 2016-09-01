perf = [49.8786 38.5571 74.7286 69.0571 64.1071 64.3714];

normalise = @(x) tertiaryop(range(x) ~= 0,(x-min(x))./(max(x)-min(x)),x);
perf = normalise(perf);

%%

objective = @(p) sum((perf-normalise(contrastGratingsMutualInformation(p(1),p(2),p(3:6)))).^2);
x0 = [1 3 1 0 0 0];
lb = [0 0 -Inf -Inf -Inf -Inf];
ub = Inf(1,6);
options = optimset(optimset(@fmincon),'Algorithm','interior-point');

params = fmincon(objective,x0,[],[],[],[],lb,ub,@contrastGratingsMutualInformationNonlinearConstraint,options);