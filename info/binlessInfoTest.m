offset = 0;
modulationDepth = logspace(0,2,5);
lambdaSs = zeros(2,8,10,5);
pss = zeros(1,8,10,5);
pids = zeros(5,10,5);

%%

preferredDirections(1,:,:) = mod(8*rand(1,10,5),8);
dPDs = mod(8*rand(1,10,5),8)-4;
preferredDirections(2,:,:) = preferredDirections(1,:,:)-dPDs;

assert(circ_rtest(pi*reshape(diff(preferredDirections),50,1)/4) > 0.05,'dPDs should be uniformly circularly distributed');
assert(circ_rtest(pi*reshape(2*abs(diff(preferredDirections)),50,1)/4) > 0.05,'abs(dPDs) should be uniformly distributed on the half circle');

%%

for gg = 1:5
    for hh = 1:10
        tic;
    %     lambdaS = randpartition([2 8],36); %[1 3 3 3 5; 1 1 3 5 5];
        PD = repmat(preferredDirections(:,hh,gg),1,8);
        x = repmat(0:7,2,1);

        lambdaS = modulationDepth(gg).*(cos(pi*(x-PD)/4)+1)/2+offset;
        lambdaSs(:,:,hh,gg) = lambdaS;

        ps = ones(1,8)/8; %randpartition(8)'; %[1 2 3 2 1]/9;
        pss(:,:,hh,gg) = ps;

        assert(abs(sum(ps)-1) < 10e-10,'p(s) must sum to 1');
        % lambda = mean(lambdaS);

        T = 1;
        deltaT = 1e-3;

        % Hr = lambda*T*(1-log(lambda));
        % Hr_s = sum(lambdaS.*T.*(1-log(lambdaS))/8,2);
        % 
        % Irs = Hr - Hr_s;
        % Irs2 = Irs/log(2);

        %%

        Hnr = zeros(2,200);
        Hlr = zeros(2,200);

        ns = size(lambdaS,2);
        Hnr_s = zeros(2,ns,200);
        Hlr_s = zeros(2,ns,200);
        Jnr_s = zeros(2,ns,200); % J => cross entropy
        Jlr_s = zeros(2,ns,200);

        for ii = 1:200
            n = ii-1;
            pn_S = exp(-lambdaS*T).*(lambdaS*T).^n/factorial(n);
            pn_S(~isfinite(pn_S)) = 0;
            pn = sum(pn_S.*repmat(ps,2,1),2);

        %     assert(max(abs(pn-mean(pn_S,2))) < 10e-10);

            if n == 0
                logku = 0;
            else
                logku = sum(log(1:n))-n*log(T);
            end

            for jj = 1:size(lambdaS,2)
                pn_s = pn_S(:,jj);
                Hnr_s(:,jj,ii) = -plogp(pn_s);
                Hlr_s(:,jj,ii) = -pn_s*logku;
                Jnr_s(:,jj,ii) = -plogq(pn_s,pn);
                Jlr_s(:,jj,ii) = -pn_s*logku;
            end

            Hnr(:,ii) = -plogp(pn);
            Hlr(:,ii) = -pn*logku;
        end

        Hlr(~isfinite(Hlr)) = 0;
        Hlr_s(~isfinite(Hlr_s)) = 0;
        Jlr_s(~isfinite(Jlr_s)) = 0;

        assert(isequal(Hlr_s,Jlr_s),'Location cross-entropy should not be stimulus-dependent');

        %%

        Hnr = sum(Hnr,2);
        Hlr = sum(Hlr,2);
        Hnr_s = sum(Hnr_s,3);
        Hnr_S = sum(Hnr_s.*repmat(ps,2,1),2);
        % assert(max(abs(Hnr_S-mean(Hnr_s,2))) < 10e-10);
        Hlr_s = sum(Hlr_s,3);
        Hlr_S = sum(Hlr_s.*repmat(ps,2,1),2);

        assert(max(abs(Hlr-Hlr_S)) < 10e-6,'Singular location entropy should not be stimulus-dependent');

        Irs = (Hnr - Hnr_S)/log(2);

        %%

        Jnr_s = sum(Jnr_s,3);
        Dkl_s = Jnr_s-Hnr_s;
        red = sum(min(Dkl_s,[],1).*ps,2)/log(2);
        % assert(abs(red-mean(min(Dkl_s,[],1),2)/log(2)) < 10e-10);

        %%

        HnR = zeros(200,200);
        HlR = zeros(200,200);
        HnR_s = zeros(200,200,ns);
        HlR_s = zeros(200,200,ns);

        for ii = 1:200
            for jj = 1:200
    %             tic;
                n = ii-1;
                m = jj-1;
                pnm_S = exp(-sum(lambdaS,1)*T).*(lambdaS(1,:)*T).^n.*(lambdaS(2,:)*T).^m/(factorial(n).*factorial(m));
                pnm_S(~isfinite(pnm_S)) = 0;

                if n == 0 && m > 0
                    logkuv_S = log(sum(1:m))-m*log(T);
                elseif n > 0 && m == 0
                    logkuv_S = log(sum(1:n))-n*log(T);
                elseif n == 0 && m == 0
                    logkuv_S = 0;
                else
                    logkuv_S = log(sum(1:n))+log(sum(1:m))-(n+m)*log(T);
                end

                HnR_s(ii,jj,:) = -plogp(pnm_S);
                HlR_s(ii,jj,:) = -pnm_S.*logkuv_S;

                pnm = sum(pnm_S.*ps,2);
        %         assert(isnan(pnm) || abs(pnm-mean(pnm_S,2)) < 10e-10);
                logkuv = logkuv_S;

                HnR(ii,jj) = -plogp(pnm);
                HlR(ii,jj) = -pnm.*logkuv;
    %             toc;
            end
        end

        HlR(~isfinite(HlR)) = 0;
        HlR_s(~isfinite(HlR_s)) = 0;

        %%

        HnR = sum(sum(HnR));
        HlR = sum(sum(HlR));
        HnR_s = squeeze(sum(sum(HnR_s,1),2));
        HnR_S = sum(HnR_s.*ps',1);
        HlR_s = squeeze(sum(sum(HlR_s,1),2));
        HlR_S = sum(HlR_s.*ps',1);

        assert(abs(HlR-HlR_S) < 10e-6,'Pair location entropy should not be stimulus-dependent');

        IRs = (HnR - HnR_S)/log(2);

        %%

        unq = Irs-red;
        syn = IRs-sum(unq)-red;

        pids(:,hh,gg) = [red;unq;syn;IRs];
        
        assert(all(isfinite(pids(:,hh,gg))),'PID must be finite');

        %%

        figure;
        bar(pids(:,hh,gg));

        toc;
        
        if any(pids(:,hh,gg) < 0.01)
            hh = hh - 1; %#ok<FXSET>
            warning('Bad PID, creating new pair!');
            continue;
        end 
    end
end

%%

% assert(all(pids(:) >= 0.01));

%%

ns = 8*ceil(logspace(1,4,7)/8);
ks = ns/8;

%%

pidh = zeros(5,10,numel(ns),10,2,5);
kld = zeros(8,2,2);
info1 = zeros(2,2);
info2 = zeros(1,2);

for gg = 1:5
    for hh = 1:10
        lambdaS = lambdaSs(:,:,hh,gg);
        ps = pss(:,:,hh,gg);

        for ii = 1:numel(ns)
            n = ns(ii);
            disp(n);
            k = ks(ii);

            for jj = 1:10
                tic;
                s = kron((1:8)',ones(k,1));
    %             s = datasample((1:8)',n,'Weights',ps);
                X = poisson_fixed_time(lambdaS(1,s)',T);
                Y = poisson_fixed_time(lambdaS(2,s)',T);

                kld(unique(s),1,1) = binlessInfo(s,X,'stratification_strategy',1,'max_embed_dim',2,'min_embed_dim',1,'kld',true);
                kld(unique(s),1,2) = binlessInfo(s,X,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1,'kld',true);

                kld(unique(s),2,1) = binlessInfo(s,Y,'stratification_strategy',1,'max_embed_dim',2,'min_embed_dim',1,'kld',true);
                kld(unique(s),2,2) = binlessInfo(s,Y,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1,'kld',true);

                info1(1,1) = binlessInfo(s,X,'stratification_strategy',1,'max_embed_dim',2,'min_embed_dim',1);
                info1(1,2) = binlessInfo(s,X,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1);

                info1(2,1) = binlessInfo(s,Y,'stratification_strategy',1,'max_embed_dim',2,'min_embed_dim',1);
                info1(2,2) = binlessInfo(s,Y,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1);

                info2(1) = binlessInfo(s,[X Y],'stratification_strategy',1,'max_embed_dim',2,'min_embed_dim',1);
                info2(2) = binlessInfo(s,[X Y],'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1);

                red = squeeze(sum(repmat(ps',[1 1 2]).*min(kld,[],2),1))';
                unq = info1-repmat(red,2,1);
                syn = info2-unq(1,:)-unq(2,:)-red;

                pidh(:,jj,ii,hh,:,gg) = [red;unq;syn;info2];
                toc;
            end
        end
    end
end

%%

save('binless_info_test.mat','lambdaSs','pss','pids','ks','info1','info2','kld','pidh','modulationDepth','preferredDirections','offset');
% assert(false);

%%

load('binless_info_test.mat'); %

%%

colours = distinguishable_colors(6);
colours = colours([1 2 3 4 6],:);

%%

apid = squeeze(median(100*abs(pidh./repmat(reshape(pids,[5 1 1 10 1 5]),[1 10 numel(ns) 1 2 1])-1),2));

%%

figure
set(gcf,'Position',[0 200 1200 800]);

for hh = 1:5
    subplot('Position',[0.05+0.33*mod(hh-1,3) 0.55-0.49*(hh > 3) 0.28 0.4]);
    hold on;

    for ii = 2

        for jj = 1:5
    %         subplot(2,3,jj);
            medianErrorbar(ks,squeeze(apid(jj,:,:,ii,hh)),2,'Color',colours(jj,:),'LineWidth',1.5);
            set(gca,'XScale','log');
            xlim([1 2500]);
        end
    end

    set(gca,'LineWidth',1.5);
    
    title(sprintf('Max Firing Rate = %3.1f Hz',modulationDepth(hh)));

    if hh == 4
        xlabel('Trials Per Stimulus');
        ylabel('Absolute Error (%)');
    end
    
    yy = ylim;
    ylim([0 yy(2)]);
end

lh = legend({'Redundancy' 'Unique 1' 'Unique 2' 'Synergy' 'Pair Information'},'Location','NorthEast','LineWidth',1.5);
set(lh,'Position',[0.05+0.33*2 0.55-0.49+0.28-0.058 0.14 0.178]);

%%

saveas(gcf,'binless_info_test_pid_abs','fig');
export_fig('binless_info_test_pid_abs','-eps','-png','-transparent','-painters');

%%

rpid = squeeze(median(100*(pidh./repmat(reshape(pids,[5 1 1 10 1 5]),[1 10 numel(ns) 1 2 1])-1),2));

%%

figure;
set(gcf,'Position',[0 200 1200 800]);

for hh = 1:5
    subplot('Position',[0.05+0.33*mod(hh-1,3) 0.55-0.49*(hh > 3) 0.28 0.4]);
    hold on;
    
    for ii = 2
        for jj = 1:5
    %         subplot(2,3,jj);
            medianErrorbar(ks,squeeze(rpid(jj,:,:,ii,hh)),2,'Color',colours(jj,:),'LineWidth',1.5);
            set(gca,'XScale','log');
            xlim([1 2500]);
        end
    end

    set(gca,'LineWidth',1.5);
    
    title(sprintf('Max Firing Rate = %3.1f Hz',modulationDepth(hh)));

    if hh == 4
        xlabel('Trials Per Stimulus');
        ylabel('Error (%)');
    end
    
%     yy = ylim;
%     ylim([0 yy(2)]);
end

lh = legend({'Redundancy' 'Unique 1' 'Unique 2' 'Synergy' 'Pair Information'},'Location','NorthEast','LineWidth',1.5);
set(lh,'Position',[0.05+0.33*2 0.55-0.49+0.28-0.058 0.14 0.178]);

%%

saveas(gcf,'binless_info_test_pid_rel','fig');
export_fig('binless_info_test_pid_rel','-eps','-png','-transparent','-painters');

%%

R = pids(1,:,:);
I1 = pids(2,:,:)+R;
I2 = pids(3,:,:)+R;
Ip = pids(5,:,:);

Is = [R;I1;I2;Ip];

%%

Rh = pidh(1,:,:,:,:,:);
I1h = pidh(2,:,:,:,:,:)+Rh;
I2h = pidh(3,:,:,:,:,:)+Rh;
Iph = pidh(5,:,:,:,:,:);

Ih = [Rh;I1h;I2h;Iph];

%%

apI = squeeze(median(100*abs(Ih./repmat(reshape(Is,[4 1 1 10 1 5]),[1 10 numel(ns) 1 2 1])-1),2));

%%

figure;
set(gcf,'Position',[0 200 1200 800]);

for hh = 1:5
    subplot('Position',[0.05+0.33*mod(hh-1,3) 0.55-0.49*(hh > 3) 0.28 0.4]);
    hold on;
    
    for ii = 2
        for jj = 1:4
    %         subplot(2,3,jj);
            medianErrorbar(ks,squeeze(apI(jj,:,:,ii,hh)),2,'Color',colours(jj,:),'LineWidth',1.5);
            set(gca,'XScale','log');
            xlim([1 2500]);
        end
    end

    set(gca,'LineWidth',1.5);
    
    title(sprintf('Max Firing Rate = %3.1f Hz',modulationDepth(hh)));

    if hh == 4
        xlabel('Trials Per Stimulus');
        ylabel('Absolute Error (%)');
    end
    
    yy = ylim;
    ylim([0 yy(2)]);
end

lh = legend({'Redundancy' 'Neuron 1 Information' 'Neuron 2 Information' 'Pair Information'},'Location','NorthEast','LineWidth',1.5);
set(lh,'Position',[0.05+0.33*2 0.55-0.49+0.28-0.022 0.175 0.142]);

%%

saveas(gcf,'binless_info_test_Ih_abs','fig');
export_fig('binless_info_test_Ih_abs','-eps','-png','-transparent','-painters');

%%

rpI = squeeze(median(100*(Ih./repmat(reshape(Is,[4 1 1 10 1 5]),[1 10 numel(ns) 1 2 1])-1),2));

%%

colours = distinguishable_colors(5);

figure;
set(gcf,'Position',[0 200 1200 800]);

for hh = 1:5
    subplot('Position',[0.05+0.33*mod(hh-1,3) 0.55-0.49*(hh > 3) 0.28 0.4]);
    hold on;
   
    for ii = 2
        for jj = 1:4
    %         subplot(2,3,jj);
            medianErrorbar(ks,squeeze(rpI(jj,:,:,ii,hh)),2,'Color',colours(jj,:),'LineWidth',1.5);
            set(gca,'XScale','log');
            xlim([1 2500]);
        end
    end

    set(gca,'LineWidth',1.5);
    
    title(sprintf('Max Firing Rate = %3.1f Hz',modulationDepth(hh)));

    if hh == 4
        xlabel('Trials Per Stimulus');
        ylabel('Error (%)');
    end
    
%     yy = ylim;
%     ylim([0 yy(2)]);
end

lh = legend({'Redundancy' 'Neuron 1 Information' 'Neuron 2 Information' 'Pair Information'},'Location','NorthEast','LineWidth',1.5);
set(lh,'Position',[0.05+0.33*2 0.55-0.49+0.28-0.022 0.175 0.142]);

%%

saveas(gcf,'binless_info_test_Ih_rel','fig');
export_fig('binless_info_test_Ih_rel','-eps','-png','-transparent','-painters');

%%

figure;
set(gcf,'Position',[0 200 1200 800]);

for hh = 1:5
    subplot('Position',[0.05+0.33*mod(hh-1,3) 0.55-0.49*(hh > 3) 0.28 0.4]);
    bar(pids(:,:,hh),'LineWidth',1.5);
    set(gca,'LineWidth',1.5,'XTickLabel',{'Red' 'Unq 1' 'Unq 2' 'Syn' 'Pair Info'});
    
    title(sprintf('Max Firing Rate = %3.1f Hz',modulationDepth(hh)));
    
    if mod(hh-1,3) == 0
        ylabel('Information (bits)');
    end
    
    ylim([0 3]);
end

R = reshape(R,50,1);
Ip = reshape(Ip,50,1);
S = reshape(pids(4,:,:),50,1);

% subplot('Position',[0.05+0.33*2 0.55-0.49 0.28 0.4]);
% box on;
% hold on;
% plot(Ip,100*R./Ip,Ip,100*S./Ip,'LineStyle','none','Marker','o');
% 
% %%
% 
% [rhoRI,pRI] = corr(Ip,R./Ip,'type','Spearman');
% [rhoSI,pSI] = corr(Ip,S./Ip,'type','Spearman');
% 
% dPDs = reshape(diff(preferredDirections,1),50,1);
% 
% %%
% 
% [rhoRD,pRD] = corr(abs(dPDs(:)),R./Ip,'type','Spearman');
% [rhoSD,pSD] = corr(abs(dPDs(:)),S./Ip,'type','Spearman');
% 
% %%
% 
% qRI = polyfit(Ip,100*R./Ip,2);
% qSI = polyfit(Ip,100*S./Ip,2);
% 
% %%
% 
% fplot(@(x) polyval(qRI,x),[0 3],'k--');
% fplot(@(x) polyval(qSI,x),[0 3],'k:');
% 
% set(gca,'LineWidth',1.5);
% set(get(gca,'Children'),'LineWidth',1.5);
% 
% legend({'Redundancy (data)' 'Synergy (data)' 'Redundancy (fit)' 'Synergy (fit)'},'LineWidth',1.5,'Location','NorthWest')
% 
% xlabel('Pair Information (bits)');
% ylabel('% of Pair Information');
% ylim([0 100]);

%%

saveas(gcf,'binless_info_test_target_pids_abs','fig');
export_fig('binless_info_test_target_pids_abs','-eps','-png','-transparent','-painters');

% assert(false);

%%

figure;
set(gcf,'Position',[100 100 1000 500]);

%%

subplot('Position',[0.06 0.1 0.425 0.85]);
hold on;
plot(Ip,100*R./Ip,Ip,100*S./Ip,'LineStyle','none','Marker','o');

%%

[bRI,~,~,~,statsRI] = regress(100*R./Ip,[ones(50,1) Ip],2);
[bSI,~,~,~,statsSI] = regress(100*S./Ip,[ones(50,1) Ip],2);

%%

fplot(@(x) bRI(1)+bRI(2)*x,[0 4],'k--');
fplot(@(x) bSI(1)+bSI(2)*x,[0 4],'k:');

set(gca,'LineWidth',1.5);
set(get(gca,'Children'),'LineWidth',1.5);

% legend({'Redundancy (data)' 'Synergy (data)' 'Redundancy (fit)' 'Synergy (fit)'},'LineWidth',1.5,'Location','NorthEast')

xlabel('Pair Information (bits)');
ylabel('% of Pair Information');
xlim([0 3]);
ylim([0 100]);

%%

subplot('Position',[0.55 0.1 0.425 0.85]);
hold on;
plot(abs(dPDs(:)),100*R./Ip,abs(dPDs(:)),100*S./Ip,'LineStyle','none','Marker','o');

%%

[bRD,~,~,~,statsRD] = regress(100*R./Ip,[ones(50,1) abs(dPDs(:))],2);
[bSD,~,~,~,statsSD] = regress(100*S./Ip,[ones(50,1) abs(dPDs(:))],2);

%%

fplot(@(x) bRD(1)+bRD(2)*x,[0 4],'k--');
fplot(@(x) bSD(1)+bSD(2)*x,[0 4],'k:');

set(gca,'LineWidth',1.5);
set(get(gca,'Children'),'LineWidth',1.5);

legend({'Redundancy (data)' 'Synergy (data)' 'Redundancy (fit)' 'Synergy (fit)'},'LineWidth',1.5,'Location','NorthEast')

xlabel('Absolute Difference in Preferred Direction');
% ylabel('% of Pair Information');
xlim([0 4]);
ylim([0 100]);

%%

saveas(gcf,'binless_info_test_pids_versus_info_dpd','fig');
export_fig('binless_info_test_pids_versus_info_dpd','-eps','-png','-transparent','-painters');

%%

yR = R./Ip;
yS = S./Ip;

X = [Ip abs(dPDs(:))];

%%

mdlR = LinearModel.fit(X,yR,'interactions','VarNames',{'Info' 'dPD' 'Red'});

%%

mdlS = LinearModel.fit(X,yS,'interactions','VarNames',{'Info' 'dPD' 'Syn'});

%%

[Imesh,Dmesh] = meshgrid(0:0.01:3,0:0.1:4);
Rpred = predict(mdlR,[Imesh(:) Dmesh(:)]);
figure;
surf(Imesh,Dmesh,reshape(Rpred,size(Imesh)));
shading interp;
zlim([0 1]);

%%

Inew = kron(mean(reshape(Ip,10,5))',ones(41,1));
Dnew = repmat((0:0.1:4)',5,1);
Rnew = predict(mdlR,[Inew Dnew]);
Snew = predict(mdlS,[Inew Dnew]);

%%

figure;
bar(100*pids(1:4,:)./repmat(pids(5,:),4,1),'LineWidth',1.5);
set(gca,'LineWidth',1.5,'XTickLabel',{'Redundancy' 'Unique 1' 'Unique 2' 'Synergy' 'Pair Information'});
ylabel('Percentage of Pair Information');

%%

saveas(gcf,'binless_info_test_target_pids_rel','fig');
export_fig('binless_info_test_target_pids_rel','-eps','-png','-transparent','-painters');