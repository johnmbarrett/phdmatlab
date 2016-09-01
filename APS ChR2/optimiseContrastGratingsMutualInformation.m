perf = [46.4786 45.7071 60.3357 53.1571 53.6500 55.2357; 49.8786 38.5571 74.7286 69.0571 64.1071 64.3714]';

normalise = @(x) tertiaryop(range(x(:)) ~= 0,(x-min(x(:)))./(max(x(:))-min(x(:))),x);

perf = normalise(perf); %max(0,perf/50-1); 

%%

% 3rd degree polynomial version
% objective = @(p) sum((perf-normalise(contrastGratingsMutualInformation(p(1),p(2),p(3:6)))).^2);
% x0 = [1 3 1 0 0 0];
% lb = [0 0 0 -Inf -Inf -Inf];
% ub = [Inf 6 Inf Inf Inf Inf];
% options = optimset(optimset(@fmincon),'Algorithm','interior-point');
% 
% params = fmincon(objective,x0,[],[],[],[],lb,ub,@contrastGratingsMutualInformationNonlinearConstraint,options);

% m free in [0,6]
% 
% objective = @(p) sum(sum((perf-normalise([                  ...
%     contrastGratingsMutualInformation(p(1),p(2),p(3:5))     ...
%     contrastGratingsMutualInformation(p(6),p(7),p(8:10))    ...
%     ])).^2));
% 
% x0 = [1 3 1 3 1 1 3 1 3 1]';
% lb = [0 0 0 0 0 0 0 0 0 0]';
% ub = [250 6 250 6 Inf 250 6 250 6 Inf]';
% A = [1 0 1 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 1 0 0];
% b = [250; 250];

% m fixed at 1

objective = @(p) sum(sum((perf(:,2)-normalise(contrastGratingsMutualInformation(p(1),p(2:4)))).^2));
    
    %     contrastGratingsMutualInformation(p(6),p(6:8))  ...

x0 = ones(4,1);
lb = zeros(4,1);
ub = [250 250 Inf Inf]; % 250 250 Inf Inf]';
A = [1 1 0 0]; % 0 0 0 0; 0 0 0 0 1 1 0 0];
b = 250; % 250];
options = optimset(optimset(@fmincon),'Algorithm','interior-point');

% params = fmincon(objective,x0,A,b,[],[],lb,ub,[],options);

problem = createOptimProblem('fmincon','Aineq',A,'bineq',b,'lb',lb,'objective',objective,'options',options,'ub',ub,'x0',x0);
gs = GlobalSearch;

params = run(gs,problem);