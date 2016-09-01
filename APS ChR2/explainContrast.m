maskResp = 0.5;
gratingResp = [0.5 0.5 0.5 0.5;     ...
               0.6 0.55 0.45 0.45;  ...
               0.85 0.75 0.5 0.5;   ...
               0.8 0.7 0.45 0.45;   ...
               0.75 0.65 0.4 0.4;   ...
               0.7 0.6 0.35 0.35]';
           
%%

mask = 160*ones(664,1);
contrasts = 0.1:0.1:0.6;

%%

gratings = zeros(664,4,6);
x = kron(1-2*mod(0:5,2)',ones(160,1));

for ii = 1:6
    for jj = 1:4
        gratings(:,jj,ii) = 160*(1+contrasts(ii)*x((jj-1)*80+(1:664)));
    end
    
%     subplot(2,3,ii);
%     plot(gratings(:,:,ii));
%     ylim([0 2]);
%     xlim([1 664]);
end

gratings = squeeze(mat2cell(gratings,664,ones(4,1),ones(6,1)));

%%

resps = [maskResp; gratingResp(:)];
images = [ones(664,1); gratings(:)];

%%

objFun = @(p) mean((resps-cellfun(@(I) simpleChR2RGCLNModel(I,1,p(1),60,p(2),p(3),p(4),p(5)),images)).^2);
% objFun = @(p,x) cellfun(@(I) simpleChR2RGCLNModel(I,1,p(1),60,p(2),p(3),p(4),p(5)),x);

%%

% x a b t g s
x0 = [332 1 0 128 100 0]';
lb = [0 0 0 0 10 0]';
ub = [664 250 250 255 1000 250]';
A = [0 1 1 0 0 1];
b = 250;

%%

% options = optimset(optimset(@lsqcurvefit),'Algorithm','interior-point');
% p = lsqcurvefit(objFun,x0,images,resps,lb,ub);

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'lb',lb,'ub',ub,'x0',x0,'options',options);
gs = GlobalSearch;
p = run(gs,problem);

%%

perf = [46.4786 45.7071 60.3357 53.1571 53.6500 55.2357; 49.8786 38.5571 74.7286 69.0571 64.1071 64.3714]';

% normalise = @(x) tertiaryop(range(x(:)) ~= 0,(x-min(x(:)))./(max(x(:))-min(x(:))),x);

perf = max(0,perf/50-1); % normalise(perf);

%%

objFun = @(p) sum((perf(:,2)-mean(contrastGratingsMutualInformationSingleNeuron(p(6),simpleChR2RGCLNModel(mask,1,p(1),30,p(2),p(3),p(4),p(5)),cellfun(@(I) simpleChR2RGCLNModel(I,1,p(1),60,p(2),p(3),p(4),p(5)),gratings)),1)').^2);

%%

% options = optimset(optimset(@lsqcurvefit),'Algorithm','interior-point');
% p = lsqcurvefit(objFun,x0,[],perf(:,2),lb,ub);

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'Aineq',A,'bineq',b,'lb',lb,'ub',ub,'x0',x0,'options',options);
gs = GlobalSearch;
p = run(gs,problem)';

%%

% s x p 
x0 = [0 332 0 0 1 0]';
lb = [0 0 0 -Inf -Inf -Inf]';
ub = [250 664 Inf Inf Inf Inf]';
nonlcon = @polyChR2RGCLNModelNonlinearConstrains;

%%

objFun = @(p) sum((perf(:,2)-mean(contrastGratingsMutualInformationSingleNeuron(p(1),polyChR2RGCLNModel(mask,1,p(2),30,p(3:6)),cellfun(@(I) polyChR2RGCLNModel(I,1,p(2),60,p(3:6)),gratings)),1)').^2);

%%

% options = optimset(optimset(@lsqcurvefit),'Algorithm','interior-point');
% p = lsqcurvefit(objFun,x0,[],perf(:,2),lb,ub);

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'nonlcon',nonlcon,'lb',lb,'ub',ub,'x0',x0,'options',options);
gs = GlobalSearch;
p = run(gs,problem)';

%%

load('contrast_example_data.mat');

%%

resps = [exampleMaskResp; reshape(exampleGratResp,1200,1)];
images = [mask; gratings(:)];

%%

% x a b t g s
x0 = [332 1 0 128 1]'; % 0]';
lb = [0 1 0 0 1]'; % 0]';
ub = [664 250 250 224 10]'; %250]';
% A = [0 255 1 -1]'; % 0 1];
% b = 250;

%%

objFun = @(p) mean((resps-[repmat(rectifyingChR2RGCLNModel(ones(664,1),1,p(1),30,p(2),p(3),p(4),p(5)),1200,1);reshape(repmat(reshape(cellfun(@(I) rectifyingChR2RGCLNModel(I,1,p(1),30,p(2),p(3),p(4),p(5)),gratings),[1 4 6]),[50 1 1]),1200,1)]).^2);

%%

% options = optimset(optimset(@lsqcurvefit),'Algorithm','interior-point');
% p = lsqcurvefit(objFun,x0,images,resps,lb,ub);

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'lb',lb,'ub',ub,'x0',x0,'options',options);
gs = GlobalSearch;
p = run(gs,problem);

%%

y = resps;
X = [repmat([2.5 0],1200,1); repmat(kron((1:4)',ones(50,1)),6,1) kron(contrasts',ones(200,1))];

%%

mdl = GeneralizedLinearModel.fit(X,y,'interactions','Distribution','poisson','VarNames',{'Phase' 'Contrast' 'Spike Count'});

%%

mu = feval(mdl,[2.5 0; repmat((1:4)',6,1) kron(contrasts',ones(4,1))]);

%%

objFun = @(p) mean((mu-cellfun(@(I) generalChR2RGCLNModel(I,1,p(1),p(2),p(3:6),@(x,q) q(1)*log(1+exp(x-q(3))).^q(4)+q(2)),images)).^2);

%%

% x s a b mu sigma
x0 = [332 30 1 0 128 1];
lb = [0 1 0 0 0 1];
ub = [664 664 250 250 255 10];
% A = [0 0 1 1 0 0; -1 1 0 0 0 0; 1 1 0 0 0 0];
A = [-1 1 0 0 0 0; 1 1 0 0 0 0];
% b = [250; 0; 664];
b = [0; 664];

%%

% options = optimset(optimset(@lsqcurvefit),'Algorithm','interior-point');
% p = lsqcurvefit(objFun,x0,images,resps,lb,ub);

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'Aineq',A,'bineq',b,'lb',lb,'ub',ub,'x0',x0,'options',options);
gs = GlobalSearch;
p = run(gs,problem);

%%

M = mask';
G = zeros(4,6,664);

for ii = 1:4
    for jj = 1:6
        G(ii,jj,:) = gratings{ii,jj};
    end
end

G = reshape(G,24,664);
        
I = [repmat(M,1200,1);kron(G,ones(50,1))];
% I = I(:,1:80:end);

%%

mdl = GeneralizedLinearModel.fit(I,y,'linear','Distribution','poisson');

%%

kernel = @(a,b,mu,sigma,x) a*((x/128-1)*normpdf((1:664)',mu,sigma)-b);
loglik = @(a,b,mu,sigma,n,x) n.*kernel(a,b,mu,sigma,x)-exp(kernel(a,b,mu,sigma,x))-factorial(n);
objFun = @(p) -sum(loglik(p(1),p(2),p(3),p(4),resps,I));

%%

x0 = [1 0 332 30];
lb = [0 -Inf 0 1];
ub = [Inf Inf 664 664];
A = [0 0 -1 2; 0 0 1 2];
b = [0; 664];

%%

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'Aineq',A,'bineq',b,'lb',lb,'ub',ub,'x0',x0,'options',options);
gs = GlobalSearch;
p = run(gs,problem);

%%

light = [ones(1200,1); kron(1+contrasts',ones(200,1))];
dark = [ones(1200,1); kron(1-contrasts',ones(200,1))];

lambda = @(a,b,t,g,x,y) a.*exp(g*(x-t))+b.*exp(g*(y-t));
loglik = @(a,b,t,g,x,y,n) n.*log(lambda(a,b,t,g,x,y))-lambda(a,b,t,g,x,y);

objFun = @(p) -sum(loglik([0.5*ones(1200,1);repmat(kron(p(1:4),ones(50,1)),6,1)],[0.5*ones(1200,1);repmat(kron(p(5:8),ones(50,1)),6,1)],p(9),p(10),light,dark,resps));

%%

%%

x0 = [1 1 1 1 0 0 0 0 0 1]';
lb = [0 0 0 0 0 0 0 0 -Inf 0]';
ub = [1 1 1 1 1 1 1 1 Inf Inf]';
A = [1 0 0 0 1 0 0 0 0 0; 0 1 0 0 0 1 0 0 0 0; 0 0 1 0 0 0 1 0 0 0; 0 0 0 1 0 0 0 1 0 0];
b = ones(4,1);

%%

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'Aeq',A,'beq',b,'lb',lb,'ub',ub,'x0',x0,'options',options);
gs = GlobalSearch;
p = run(gs,problem);

%%

bwgratings = cellfun(@(I) repmat((I-min(I(:)))./(max(I(:))-min(I(:))),1,664),gratings(:,end),'UniformOutput',false);

%%

[Y,X] = ndgrid(1:664,1:664);

%%

isInEllipse = @(x,y,a,b,theta) ((X-x)*cos(theta)+(Y-y)*sin(theta)).^2/a^2 + ((Y-y)*cos(theta)-(X-x)*sin(theta)).^2/b^2 <= 1;

%%

objFun = @(q) mean((p(1:4)-cellfun(@(I) sum(sum(logical(I) & isInEllipse(q(1),q(2),q(3),q(4),q(5))))/(664*644),bwgratings)).^2);

%%

y0 = [332 332 30 30 0]';
lb = [0 0 1 1 0]';
ub = [664 664 332 332 2*pi]';

%%

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'lb',lb,'ub',ub,'x0',y0,'options',options);
gs = GlobalSearch;
q = run(gs,problem);

%%

objFun = @(q) mean((p(1:4)-cellfun(@(I) sum(sum(I.*gauss2d(X-q(1),Y-q(2),q(3),q(4),q(5),2))),bwgratings)).^2);

%%

y0 = [332 332 30 30 0]';
lb = [0 0 1 1 0]';
ub = [664 664 332 332 2*pi]';

%%

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'lb',lb,'ub',ub,'x0',y0,'options',options);
gs = GlobalSearch;
q = run(gs,problem);

%%

objFun = @(q) mean((p(1:4)-cellfun(@(I) sum(I(:,1).*gauss1d((1:664)'-q(1),q(2),2)),bwgratings)).^2);

%%

y0 = [332 30]';
lb = [0 1]';
ub = [664 332]';
A = [-1 2; 1 2];
b = [0; 664];

%%

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'Aineq',A,'bineq',b,'lb',lb,'ub',ub,'x0',y0,'options',options);
gs = GlobalSearch;
q = run(gs,problem);

%%

% lambda = @(a,b,t,g,x,y) a.*exp(g*(x-t))+b.*exp(g*(y-t));
I = contrastGratingsMutualInformationSingleNeuron(0,lambda(1,0,p(9),p(10),1,1),lambda(  ...
    repmat(p(1:4),6,1),             ... a
    repmat(p(5:8),6,1),             ... b
    p(9), p(10),                    ... t,g
    kron(1+contrasts',ones(4,1)),   ... x
    kron(1-contrasts',ones(4,1))    ... y
    ));

%%

rectify = @(a,b,x) a*max(0,x-b);
objFun = @(r) mean((perf(:,2)-contrastGratingsMutualInformationSingleNeuron(r(1),rectify(r(2),r(3),1),rectify(r(2),r(3),r(4)*(1+contrasts')+(1-r(4))*(1-contrasts')))).^2);

%%

z0 = [0 1 0 0.5];
lb = [0 0 0 0];
ub = [250 250 2 1];

%%

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'lb',lb,'ub',ub,'x0',z0,'options',options);
gs = GlobalSearch;
r = run(gs,problem);

%%

objFun = @(r) mean((perf(:,2)-contrastGratingsMutualInformationSingleNeuron(r(1),rectify(r(2),r(3),1),r(4)*rectify(r(2),r(3),1+contrasts')+(1-r(4))*rectify(r(2),r(3),1-contrasts'))).^2);

%%

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'lb',lb,'ub',ub,'x0',z0,'options',options);
gs = GlobalSearch;
r = run(gs,problem);

%%

objFun = @(r) mean((perf(:,2)-contrastGratingsMutualInformation(r(1),[],[],[rectify(r(2),r(3),1) rectify(r(5),r(6),1)],[r(4)*rectify(r(2),r(3),1+contrasts')+(1-r(4))*rectify(r(2),r(3),1-contrasts') r(7)*rectify(r(5),r(6),1+contrasts')+(1-r(7))*rectify(r(5),r(6),1-contrasts')],1856)).^2);

%%

z0 = [0 1 0 0.5 1 0 0.5];
lb = [0 0 0 0 0 0 0];
ub = [250 250 2 1 250 2 1];

%%

options = optimset(optimset(@fmincon),'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',objFun,'lb',lb,'ub',ub,'x0',z0,'options',options);
gs = GlobalSearch;
r = run(gs,problem);

%%

sigma = r(1);
mu = [rectify(r(2),r(3),1) rectify(r(5),r(6),1)];
lambda = [r(4)*rectify(r(2),r(3),1+contrasts')+(1-r(4))*rectify(r(2),r(3),1-contrasts') r(7)*rectify(r(5),r(6),1+contrasts')+(1-r(7))*rectify(r(5),r(6),1-contrasts')];

%%

I1 = contrastGratingsMutualInformationSingleNeuron(sigma,mu(1),lambda(:,1));
I2 = contrastGratingsMutualInformationSingleNeuron(sigma,mu(2),lambda(:,2));
Ip = contrastGratingsMutualInformation(sigma,[],[],mu,lambda,1856);

%%

figure;
hs = plot(100*contrasts,[I1 I2 Ip],'LineWidth',1.5);
box off;
set(gca,'LineWidth',1.5,'XTick',10:10:60,'YTick',0:0.1:0.6);
xlabel('Michelson Contrast (%)');
xlim([10 60]);
ylabel('Information (bits)');
ylim([0 0.6]);

pos = get(gca,'Position');
colorOrder = get(gca,'ColorOrder');

ax = axes;
hs(end+1) = plot(100*contrasts,50+50*perf(:,2),'Color',colorOrder(4,:),'LineWidth',1.5);
box off;
set(ax,'Color','none','LineWidth',1.5,'Position',pos,'XAxisLocation','top','XTick',[],'YAxisLocation','right','YColor',colorOrder(4,:));
xlim([10 60]);
yh = ylabel('Decoder Performance (%)');
set(yh,'Rotation',270,'VerticalAlignment','bottom');
ylim([50 80]);

legend(hs([1 2 3 4]),{'Model Neuron 1' 'Model Neuron 2' 'Model Neuron Pair' 'Decoder Performance'},'LineWidth',1.5,'Location','NorthEast');

%%

saveas(gcf,'new_contrast_model_info_vs_perf','fig');
export_fig('new_contrast_model_info_vs_perf','-eps','-png',-'m4','-transparent','-painters');

%%

figure;
hold on;
h = ezplot(@(x) 4*(rectify(r(2),r(3),x)+r(1)),[0 2]);
set(h,'Color',colorOrder(1,:));
h = ezplot(@(x) 4*(rectify(r(5),r(6),x)+r(1)),[0 2]);
set(h,'Color',colorOrder(2,:));
line([1 1],ylim,'Color','k','LineStyle','--');
set(gca,'LineWidth',1.5,'XTick',0:0.2:2,'XTickLabel',0:20:200);
set(get(gca,'Children'),'LineWidth',1.5);
title('');
xlabel('Irradiance (% of Mask)');
ylabel('Mean Firing Rate (Hz)');
legend({'Neuron 1' 'Neuron 2' 'Mask Irradiance'},'LineWidth',1.5,'Location','NorthWest');

%%

saveas(gcf,'new_contrast_model_responses','fig');
export_fig('new_contrast_model_responses','-eps','-png',-'m4','-transparent','-painters');

%%

figure;
set(gcf,'Position',[32 300 1000 500]);
yy = [-0.5 4.5; 27 60];

hs = zeros(3,1);

for ii = 1:2
    subplot('Position',[0.06+0.49*(ii-1) 0.09 0.44 0.85]);
    hold on;
    
    m = poissinv(0.50,lambda(:,ii));
    l = m-poissinv(0.25,lambda(:,ii));
    u = poissinv(0.75,lambda(:,ii))-m;
    
    hs(2) = line([0 70],[1 1]*poissinv(0.5,mu(ii)+sigma),'Color',colorOrder(2,:),'LineStyle','-','LineWidth',1.5);
    hs(3:4) = line([0 0; 70 70],[1;1]*poissinv([0.25 0.75],mu(ii)+sigma),'Color',colorOrder(2,:),'LineStyle','--','LineWidth',1.5);
    hs(5:6) = line([0 0; 70 70],[1;1]*poissinv([0.025 0.975],mu(ii)+sigma),'Color',colorOrder(2,:),'LineStyle',':','LineWidth',1.5);
    
    hs(1) = errorbar(100*contrasts',m,l,u,'LineWidth',1.5);
    
    set(gca,'LineWidth',1.5);
    
    title(sprintf('Neuron %d',ii));
    
    xlabel('Michelson Contrast (%)');
    xlim([0 70]);
    
    if ii == 1
        legend(hs([1:3 5]),{'Grating' 'Mask (median)' 'Mask (IQR)' 'Mask (95% CI)'},'LineWidth',1.5,'Location','NorthWest');
        ylabel('# Spikes');
    end
    
    ylim(yy(ii,:));
end

%%

saveas(gcf,'new_contrast_model_grating_vs_mask','fig');
export_fig('new_contrast_model_grating_vs_mask','-eps','-png',-'m4','-transparent','-painters');

%%

figure;
fplot(@(x) contrastGratingsMutualInformation(r(1),[],[],mu,[r(4)*rectify(r(2),r(3),1+x(:))+(1-r(4))*rectify(r(2),r(3),1-x(:)) r(7)*rectify(r(5),r(6),1+x(:))+(1-r(7))*rectify(r(5),r(6),1-x(:))]),[0 1]);

%%

fs = get(gca,'Children');
set(fs,'Color',colorOrder(3,:),'LineWidth',1.5);

%%

hold on;
fs(2) = ezplot(@(x) contrastGratingsMutualInformationSingleNeuron(r(1),mu(1),r(4)*rectify(r(2),r(3),1+x(:))+(1-r(4))*rectify(r(2),r(3),1-x(:))),[0 1]);
fs(3) = ezplot(@(x) contrastGratingsMutualInformationSingleNeuron(r(1),mu(2),r(7)*rectify(r(5),r(6),1+x(:))+(1-r(7))*rectify(r(5),r(6),1-x(:))),[0 1]);

%%

set(fs(2),'Color',colorOrder(1,:),'LineWidth',1.5);
set(fs(3),'Color',colorOrder(2,:),'LineWidth',1.5);

%%

set(gca,'LineWidth',1.5,'XTick',0:0.1:1,'XTickLabel',0:10:100);
title('');
xlabel('Michelson Contrast (%)');
xlim([0.1 1]);
ylabel('Information (bits)');
ylim([0 0.65]);
legend(fs([2 3 1]),{'Model Neuron 1' 'Model Neuron 2' 'Model Neuron Pair'},'LineWidth',1.5,'Location','NorthWest');

%%

saveas(gcf,'new_contrast_model_info_vs_contrast','fig');
export_fig('new_contrast_model_info_vs_contrast','-eps','-png',-'m4','-transparent','-painters');