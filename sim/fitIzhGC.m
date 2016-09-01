function [params,fval,exitflag,output,solutions] = fitIzhGC
    % estimated using GlobalSearch with fmincon (see below)
    a = 0.0115;
    b = 0.0584;
    c = -56.0346;
    d = 6.7354;

    tI = [0 499 500 2499 2500 2999]';
    I = [0 0 1 1 0 0]';
    amps = 15:5:50;
    targetFR = 2*0.6*amps-10;
    
    function fr = modelFn(params,amps)
        s = izhNeuron(-60,2999,tI,repmat(I,1,numel(amps)).*repmat(amps,numel(I),1),params(1),params(2),params(3),params(4));
        
        nI = size(s,2);
        firingStart = zeros(1,nI);
        firingEnd = zeros(1,nI);
        
        for ii = 1:nI
            if sum(s(:,ii)) == 0
                firingStart(ii) = 0;
                firingEnd(ii) = 0;
                continue;
            end
            
            firingStart(ii) = find(s(:,ii),1)/1000;
            firingEnd(ii) = find(s(:,ii),1,'last')/1000;
        end
        
        responseLength = firingEnd-firingStart;
        
        fr = sum(s)./responseLength;
        
        fr(~isfinite(fr)) = 0;
    end

%     params = nlinfit(amps,targetFR,@modelFn,[a b c d]);
%     params = lsqcurvefit(@modelFn,[0.2 0.02 -65 8],amps,targetFR);
    lb = [0 0 -70 0];
    ub = [0.2 0.5 0 20];
%     options = optimset('Algorithm','active-set');
    problem = createOptimProblem('lsqcurvefit','lb',lb,'ub',ub,'objective',@modelFn,'x0',[0.2 0.02 -65 8],'xdata',amps,'ydata',targetFR); %,'options',options);
    ms = MultiStart;
    [params,fval,exitflag,output,solutions] = run(ms,problem,100);
end

%     tI = [0 499 500 2499 2500 2999]';
%     I = [0 0 1 1 0 0]';
%     amps = 15:5:50;
%     I = repmat(I,1,numel(amps)).*repmat(amps,numel(I),1);
%     targetFR = 2*0.6*amps-10;
%     
%     params = [0.2 0.02 -65 8];
%     lb = [0 0 -70 0];
%     ub = [0.2 0.5 0 20];
%     options = optimset('Algorithm','active-set');
%     
%     function mse = modelFn(params)
%         s = izhNeuron(-60,2999,tI,I,params(1),params(2),params(3),params(4));
%         
%         nI = size(s,2);
%         firingStart = zeros(1,nI);
%         firingEnd = zeros(1,nI);
%         
%         for ii = 1:nI
%             if sum(s(:,ii)) == 0
%                 firingStart(ii) = 0;
%                 firingEnd(ii) = 0;
%                 continue;
%             end
%             
%             firingStart(ii) = find(s(:,ii),1)/1000;
%             firingEnd(ii) = find(s(:,ii),1,'last')/1000;
%         end
%         
%         responseLength = firingEnd-firingStart;
%         
%         fr = sum(s)./responseLength;
%         
%         fr(~isfinite(fr)) = 0;
%         
%         mse = sum((targetFR-fr).^2);
%     end
% 
% %     params = fmincon(@modelFn,params,[],[],[],[],lb,ub,[],options);
%     problem = createOptimProblem('fmincon','lb',lb,'ub',ub,'objective',@modelFn,'x0',params,'options',options);
%     gs = GlobalSearch;
%     params = run(gs,problem);
% end

% clear all;
% load('izhdata.mat');
% bw = 50;
% kernel = normpdf(-5*bw:bw*5,0,bw);
% s = conv(s,kernel,'same');
% params = [0.2 0.02 -65 8];
% lb = [0 0 -70 0];
% ub = [0.2 0.5 0 20];
% tmax = numel(s)-1;
% tI(end) = tmax;
% I = 1.25*(I-5);
% I = interp1(tI,I,0:tmax);
% tI = 0:tmax;
% modelFn = @(params) mean((s-conv(izhNeuron(-65,tmax,tI,I,params(1),params(2),params(3),params(4)),kernel,'same')).^2);
% options = optimset('Algorithm','interior-point');
% params = fmincon(modelFn,params,[],[],[],[],lb,ub,[],options);
% [t,v,u] = izhNeuron(-65,tmax,tI,I,params(1),params(2),params(3),params(4));
% %%
% figure;
% plotyy(tI,v,tI,I);
% ylim([min(v) max(v)]);
% %%
% figure;
% t = conv(t,kernel,'same');
% plot(tI,[s t]);