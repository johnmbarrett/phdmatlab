function rd1model
% function fr = rd1model
% function [params,v,u,s,tINL,vINL,Isyn] = rd1model(targetFR)
    [tINL,yINL] = morrisLecarNetwork([0 1000],1);
    vINL = yINL(:,1:3);
    
%     modelFn = @(params) abs(targetFR-simulateGC(params,tINL,vINL));
%     
%     options = optimset('Algorithm','active-set');
%     params = fmincon( ...
%         modelFn,                        ...
%         [-30 -30 -30 10 10 10 1 1 20],  ...
%         [], [], [], [],                 ...
%         [-Inf(1,3) zeros(1,6)],         ... 
%         Inf(1,9),                       ...
%         [],                             ...
%         options                         ...
%         );
    
%     fr = zeros(50,50);
%     
%     for ii = 1:50
%         for jj = 1:50
%             tic;
            params = [-30 -30 -30 10 10 10 10 10 50];
            Isyn = getSynapticCurrent(params,vINL);
            v = izhNeuron(-65,1000,tINL,sum(Isyn,2),0.0115,0.0584,-56.0346,6.7354);
%             fr(ii,jj) = sum(s);
%             toc;
%         end
%     end

    figure
    plot(tINL,Isyn);
    hold on;
    plot(tINL,sum(Isyn,2),'Color','m');
    
    figure;
    plot(tINL,vINL);
    hold on;
    plot(v,'Color','m');
end

function fr = simulateGC(params,tINL,vINL)
    Isyn = getSynapticCurrent(params,vINL);
    
    [~,~,s] = izhNeuron(-65,1000,tINL,sum(Isyn,2));
    
    fr = sum(s);
end

function Isyn = getSynapticCurrent(params,vINL)
    n = size(vINL,1);
    threshold = repmat(params(1:3),n,1);
    gain = repmat(params(4:6),n,1);
    maxResponse = repmat(params(7:9),n,1);
    
    Isyn = maxResponse./(1+exp(-(vINL-threshold)./gain));
    Isyn = Isyn.*repmat([-1 -1 1],n,1);
end