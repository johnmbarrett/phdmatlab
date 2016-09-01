s = kron((1:8)',ones(125,1));

pdfnames = {'Beta shape' 'Beta scale' 'Chi-squared' 'Exponential' 'EV location' 'EV scale' 'Gamma shape' 'Gamma scale' 'Lognormal location' 'Lognormal scale' 'Normal location' 'Normal scale'};
pdffuns = { ...
    @(x,a) betapdf(x,a,1)           ...
    @(x,b) betapdf(x,1,b)           ...
    @chi2pdf                        ... 
    @exppdf                         ...
    @(x,mu) evpdf(x,mu,1)           ...
    @(x,sigma) evpdf(x,0,sigma)     ...
    @(x,a) gampdf(x,a,1)            ...
    @(x,b) gampdf(x,1,b)            ...
    @(x,mu) lognpdf(x,mu,1)         ...
    @(x,sigma) lognpdf(x,0,sigma)   ...
    @(x,mu) normpdf(x,mu,1)         ...
    @(x,sigma) normpdf(x,0,sigma)   ...
    };

Hfuns = { ...
    @(a) log(beta(a,1))-(a-1).*psi(a)+(a-1).*psi(a+1)   ...
    @(b) log(beta(1,b))-(b-1).*psi(b)+(b-1).*psi(1+b)   ...
    @(k) k/2+log(2*gamma(k/2))+(1-k/2).*psi(k/2)        ...
    @(m) 1+log(m)                                       ...
    @(m) 1-psi(ones(size(m)))                           ...
    @(s) log(s)-psi(1)+1                                ...
    @(a) a+log(gamma(a))+(1-a).*psi(a)                  ...
    @(b) 1+log(b)+log(gamma(1))                         ...
    @(m) 1/2+log(2*pi)/2+m                              ...
    @(s) 1/2+log(2*pi.*s.^2)/2                          ...
    @(m) log(2*pi*exp(ones(size(m))))/2                 ...
    @(s) log(2*pi*exp(1)*s.^2)/2                        ...
    };

rndfuns = { ...
    @(a,varargin) betarnd(a,1,varargin{:})          ...
    @(b,varargin) betarnd(1,b,varargin{:})          ...
    @chi2rnd                                        ... 
    @exprnd                                         ...
    @(mu,varargin) evrnd(mu,1,varargin{:})          ...
    @(sigma,varargin) evrnd(0,sigma,varargin{:})    ...
    @(a,varargin) gamrnd(a,1,varargin{:})           ...
    @(b,varargin) gamrnd(1,b,varargin{:})           ...
    @(mu,varargin) lognrnd(mu,1,varargin{:})        ...
    @(sigma,varargin) lognrnd(0,sigma,varargin{:})  ...
    @(mu,varargin) normrnd(mu,1,varargin{:})        ...
    @(sigma,varargin) normrnd(0,sigma,varargin{:})  ...
    };
    
xmins = [0 0 0   0   -Inf -Inf 0   0   0   0   -Inf -Inf];
xmaxs = [1 1 Inf Inf  Inf  Inf Inf Inf Inf Inf  Inf  Inf];

%%
pdfs = cell(12,1);
Hy = zeros(12,1);
Hyx = zeros(12,1);
D = zeros(8,12);

for ii = 1:12
    pdffun = pdffuns{ii};
    
    pdf = @(x) (pdffun(x,1)+pdffun(x,2)+pdffun(x,3)+pdffun(x,4)+pdffun(x,5)+pdffun(x,6)+pdffun(x,7)+pdffun(x,8))/8;
    pdfs{ii} = pdf;
    
    Hy(ii) = -integral(@(x) plogp(pdf(x)),xmins(ii),xmaxs(ii));
    Hyx(ii) = mean(Hfuns{ii}(1:8));
    
    for jj = 1:8
        D(jj,ii) = integral(@(x) plogq(pdffun(x,jj),pdffun(x,jj)./pdf(x)),xmins(ii),xmaxs(ii));
    end
end

I = (Hy-Hyx)/log(2);
D = D/log(2);

%%

ks = round(logspace(2,4,11)/8);

%%

% Ihat = zeros(10,11,12);
% Dhat = zeros(8,10,11,12,2);

for ii = 1:11
    k = ks(ii);
    s = kron((1:8)',ones(k,1));
    R = zeros(size(s));
    
    for jj = 1:12
        rndfun = rndfuns{jj};
        
        for kk = 1:10
            fprintf('Distribution: %s, N = %d, iteration #%02d ...\n',pdfnames{jj},k*8,kk);
            
            tic;
            for ll = 1:8
                R(k*(ll-1)+1:k*ll) = rndfun(ll,k,1);
            end
            fprintf('Generated random sample in %f seconds\n',toc);
            
            tic;
            for ll = 1:8
                P = rndfun(ll,k,1);
                Dhat(ll,kk,ii,jj,2) = mixedKLDivergence(P,R,false,1,'both',false);
            end
%             Ihat(kk,ii,jj) = binlessInfo(s,R,'stratification_strategy',0,'max_embed_dim',2,'min_embed_dim',1,'starttime',min(0,floor(min(R))),'triallength',ceil(max(R))); %victorMutualInformation(s,R,false,true);
%             Dhat(:,kk,ii,jj,1) = binlessInfo(s,R,'stratification_strategy',0,'max_embed_dim',2,'min_embed_dim',1,'starttime',min(0,floor(min(R))),'triallength',ceil(max(R))); %victorMutualInformation(s,R,false,true);
            fprintf('Calculated mutual information in %f seconds\n\n',toc);
        end
    end
end

%%

II = repmat(reshape(I,[1 1 12]),10,11);

%%

dI = abs(Ihat-II);
mdI = squeeze(mean(dI));
sdI = squeeze(std(dI));

figure;
errorbar(repmat(8*ks',1,12),mdI,sdI/sqrt(10));
set(gca,'XScale','log','XTick',8*ks);
xlim([0 1e4]);
xlabel('N');
ylabel('Absolute Error (bits)');
legend(pdfnames);

%%

pI = 100*abs(Ihat-II)./II;
% mpI = squeeze(mean(pI));
% spI = squeeze(std(pI));

figure;
hold on;
set(gca,'ColorOrder',distinguishable_colors(12));
medianErrorbar(ks,pI,'LineWidth',1.5);
set(gca,'LineWidth',1.5,'XScale','log','XTick',ks);
xlim([0 1300]);
xlabel('Trials Per Stimulus');
ylabel('Absolute Error (%)');
legend(pdfnames,'LineWidth',1.5);

%%

saveas(gcf,'binless_info_test_continuous_abs','fig');
export_fig('binless_info_test_continuous_abs','-eps','-png','-transparent','-painters');

%%

DD = repmat(reshape(D,[8 1 1 12 1]),[1 10 11 1 2]);

%%

dD = abs(Dhat-DD);
mdD = squeeze(mean(dD));
sdD = squeeze(std(dD));

figure;
errorbar(repmat(8*ks',1,12),mdD,sdD/sqrt(10));
set(gca,'XScale','log','XTick',ks);
xlim([0 1e4]);
xlabel('N');
ylabel('Absolute Error (bits)');
legend(pdfnames);

%%

pD = 100*abs(Dhat-DD)./DD;
% mpI = squeeze(mean(pI));
% spI = squeeze(std(pI));

%%

suffixes = {'subsample' 'sepsample'};

for ii = 1:2
    figure;
    hold on;
    set(gca,'ColorOrder',distinguishable_colors(12));
    medianErrorbar(ks,squeeze(median(pD(:,:,:,:,ii),1)),'LineWidth',1.5);
    set(gca,'LineWidth',1.5,'XScale','log','XTick',ks);
    xlim([0 1300]);
    xlabel('Trials Per Stimulus');
    ylabel('Absolute Error (%)');
    legend(pdfnames,'LineWidth',1.5);

    figName = sprintf('binless_kld_test_continuous_abs_%s',suffixes{ii});
    saveas(gcf,figName,'fig');
    export_fig(figName,'-eps','-png','-transparent','-painters');
end

%%

%%

ps = rand(10,8);

Hy2 = zeros(12,10);
Hyx2 = zeros(12,10);

for ii = 1:12
    Hys = Hfuns{ii}(1:8);
    
    for jj = 1:10
        p = ps(jj,:);
        
        Hyx2(ii,jj) = mean(-plogp(p)-plogp(1-p)+(1-p).*Hys);
        
        q = mean(p);
        
        Hy2(ii,jj) = -plogp(q)-plogp(1-q)+(1-q)*Hy(ii);
    end
end

I2 = (Hy2-Hyx2)/log(2);

%%

% ks = round(logspace(2,4,11)/8);
Ihat2 = zeros(10,10,12);
k = 125;

s = kron((1:8)',ones(k,1));
    
for ii = 1:12
    rndfun = rndfuns{ii};
    
    for jj = 1:10
        p = ps(jj,:);

        for kk = 1:10
            fprintf('Cont Dist: %s, instantiation = %d, iteration #%02d ...\n',pdfnames{ii},jj,kk);

            tic;
            R = zeros(size(s));
            
            for ll = 1:8
                X = rand(k,1) > p(ll);
                m = sum(X);
                r = rndfun(ll,m,1);
                R(k*(ll-1)+find(X)) = r; %num2cell(r);
                R(k*(ll-1)+find(~X)) = Inf;
            end
            fprintf('Generated random sample in %f seconds\n',toc);

            tic;
            Ihat2(kk,jj,ii) = binlessInfo(s,R,'stratification_strategy',0,'max_embed_dim',2,'min_embed_dim',1,'starttime',min(0,floor(min(R))),'triallength',ceil(max(R))); %victorMutualInformation(s,R,false,true);
            fprintf('Calculated mutual information in %f seconds\n\n',toc);
        end
    end
end

%%

II2 = repmat(reshape(I2',[1 10 12]),[10 1 1]);

%%

dI2 = abs(Ihat2-II2);
mdI2 = squeeze(mean(mean(dI2,1),2));
sdI2 = squeeze(std(mean(dI2,1),[],2));

figure;
barwitherr(sdI2/sqrt(10),mdI2);
% set(gca,'XScale','log','XTick',8*ks);
% xlim([0 1e4]);
% xlabel('N');
ylabel('Absolute Error (bits)');
% legend(pdfnames);

%%

pI2 = 100*abs(Ihat2-II2)./II2;
mpI2 = squeeze(mean(mean(pI2,1),2));
spI2 = squeeze(std(mean(pI2,1),[],2));

figure;
barwitherr(spI2/sqrt(10),mpI2);
% set(gca,'XScale','log','XTick',8*ks);
% xlim([0 1e4]);
% xlabel('N');
ylabel('Relative Error (%)');
% legend(pdfnames);

%%

Ihat2 = zeros(10000,12);

s = kron((1:8)',ones(125,1));
R = zeros(size(s));

for ii = 1:12
    rndfun = rndfuns{ii};
    
    for jj = 1:10000
        fprintf('Distribution: %s, iteration #%05d\n',pdfnames{ii},jj);
        
        tic;
        for kk = 1:8
            R(125*(kk-1)+1:125*kk) = rndfun(kk,125,1);
        end
        fprintf('Generated random sample in %f seconds\n',toc);
        
        tic;
        Ihat2(jj,ii) = victorMutualInformation(s,R,false,true);
        fprintf('Calculated mutual information in %f seconds\n\n',toc);
    end
end 

%%

ns = ceil(logspace(0.5,3,6))';

dI = zeros(10,6,12);

for ii = 1:6
    n = ns(ii);
    idx = randperm(10000);
    
    for jj = 1:10
        sample = Ihat2(idx(n*(jj-1)+1:n*jj),:);
        sampleMean = squeeze(mean(sample))';
        
        dI(jj,ii,:) = abs(sampleMean-I);
    end
end

pI = 100*dI./repmat(reshape(I,[1 1 12]),10,6);

%%

mdI = squeeze(mean(dI));
sdI = squeeze(std(dI));
figure;
errorbar(repmat(ns,1,12),mdI,sdI);
set(gca,'XScale','log','XTick',ns);
xlim([0 1e3]);
xlabel('N');
ylabel('Absolute Error (nats)');
legend(pdfnames);

%%

mpI = squeeze(mean(pI));
spI = squeeze(std(pI));

figure;
errorbar(repmat(ns,1,12),mpI,spI);
set(gca,'XScale','log','XTick',ns);
xlim([0 1e3]);
xlabel('N');
ylabel('Relative Error (%)');
legend(pdfnames);

%%

mus = kron((1:8)',ones(1,2));
sigma = [1 0.5; 0.5 1];

pdf = @(x) (mvnpdf(x,mus(1,:),sigma) + mvnpdf(x,mus(2,:),sigma) + mvnpdf(x,mus(3,:),sigma) + mvnpdf(x,mus(4,:),sigma) + mvnpdf(x,mus(5,:),sigma) + mvnpdf(x,mus(6,:),sigma) + mvnpdf(x,mus(7,:),sigma) + mvnpdf(x,mus(8,:),sigma))/8;

Hy2 = -integral2(@(x,y) plogp(pdf([x(:) y(:)]))',-Inf,Inf,-Inf,Inf);
Hyx2 = 1+log(2*pi)+log(det(sigma))/2;
I2 = Hy2-Hyx2;

%%

Ihat3 = zeros(10,11);

for ii = 1:11
    k = ks(ii);
    s = kron((1:8)',ones(k,1));
    R = zeros(8*k,2);
    
    for jj = 1:10
        for kk = 1:8
            R(k*(kk-1)+1:k*kk,:) = mvnrnd(mus(kk,:),sigma,k);
        end
        
        Ihat3(jj,ii) = victorMutualInformation(s,R,false,true);
    end
end

%%

dI = Ihat3-I2;
mdI = squeeze(mean(dI));
sdI = squeeze(std(dI));

figure;
errorbar(8*ks,mdI,sdI/sqrt(10));
set(gca,'XScale','log','XTick',8*ks);
xlim([0 1e4]);
xlabel('N');
ylabel('Absolute Error (nats)');

%%

pI = 100*(Ihat3-I2)./I2;
mpI = squeeze(mean(pI));
spI = squeeze(std(pI));

figure;
errorbar(8*ks,mpI,spI/sqrt(10));
set(gca,'XScale','log','XTick',8*ks);
xlim([0 1e4]);
xlabel('N');
ylabel('Relative Error (%)');

%%

Ihat4 = zeros(10000,1);
s = kron((1:8)',ones(125,1));
R = zeros(1000,2);

for ii = 1:10000
    fprintf('Multivariate normal, iteration #%d\n',ii);
    
    tic;
    for jj = 1:8
        R(125*(jj-1)+1:125*jj,:) = mvnrnd(mus(jj,:),sigma,125);
    end
    fprintf('Generated random sample in %f seconds\n',toc);
    
    tic;
    Ihat4(ii) = victorMutualInformation(s,R,false,true);
    fprintf('Calculated mutual information in %f seconds\n',toc);
end

%%

dI2 = zeros(10,6);

for ii = 1:6
    n = ns(ii);
    idx = randperm(10000);
    
    for jj = 1:10
        sample = Ihat4(idx(n*(jj-1)+1:n*jj),:);
        sampleMean = squeeze(mean(sample))';
        
        dI2(jj,ii) = abs(sampleMean-I2);
    end
end

pI2 = 100*dI2/I2;

%%

mdI = squeeze(mean(dI2));
sdI = squeeze(std(dI2));
figure;
errorbar(ns,mdI,sdI);
set(gca,'XScale','log','XTick',ns);
xlim([0 1e3]);
xlabel('N');
ylabel('Absolute Error (nats)');

%%

mpI = squeeze(mean(pI2));
spI = squeeze(std(pI2));

figure;
errorbar(ns,mpI,spI);
set(gca,'XScale','log','XTick',ns);
xlim([0 1e3]);
xlabel('N');
ylabel('Relative Error (%)');

%%

s = repmat(1:8,4,1);
k = repmat((0:3)',1,8);

nCk = zeros(4,8);

for ii = 1:4
    nCk(ii,:) = nchoosek(3,ii-1);
end

pns = nCk.*(s/9).^k.*(1-s/9).^(3-k);
Hns = -sum(sum(plogp(pns)))/8;

pn = sum(pns/8,2);
Hn = -sum(plogp(pn));

In = Hn-Hns;

%%

Htns = sum(sum((pns/8).*(k/2).*(1+log(2*pi))));
Htn = zeros(4,1);
integraln = {@integral,@integral2,@integral3};

%%

for ii = 2:4
    x = ii-1;
    mu = ones(1,x);
    sigma = eye(x);
    ftn = @(x) (mvnpdf(x,mu,sigma)+mvnpdf(x,2*mu,sigma)+mvnpdf(x,3*mu,sigma)+mvnpdf(x,4*mu,sigma)+mvnpdf(x,5*mu,sigma)+mvnpdf(x,6*mu,sigma)+mvnpdf(x,7*mu,sigma)+mvnpdf(x,8*mu,sigma))/8;
    limits = num2cell((1-2*mod(1:2*x,2)).*Inf(1,2*x));
    Htn(ii) = -integraln{x}(@(varargin) plogp(ftn(vertcat(varargin{:})'))',limits{:});
end

%%
    
Htn = sum(Htn.*pn);

It = Htn-Htns;

Itn = In+It;

%%

Ihat5 = zeros(10,11);

for ii = 1:11
    k = ks(ii);
    s = kron((1:8)',ones(k,1));
    X = zeros(8*k,1);
    Y = zeros(8*k,3);
    
    for jj = 1:10
        fprintf('Mixed-pair RV, n = %d, iteration #%d...\n',8*k,jj);
        
        tic;
        for kk = 1:8
            X(k*(kk-1)+1:k*kk) = binornd(3,kk/9,k,1);
            
            for ll = 1:3
                idx = find(X(k*(kk-1)+1:k*kk) == ll);
                mu = kk*ones(1,ll);
                sigma = eye(ll);
                Y(k*(kk-1)+idx,1:ll) = mvnrnd(mu,sigma,numel(idx));
            end
        end
        
        S = cell(1,4);
        R = cell(1,4);
        
        for kk = 1:4
            idx = X == kk-1;
            S{kk} = s(idx);
            R{kk} = Y(idx,1:max(1,kk-1));
        end
        fprintf('Generated random sample in %f seconds\n',toc);
        
        tic;
        Ihat5(jj,ii) = victorMutualInformation(S,R,[1 0 0 0],true);
        fprintf('Calculated mutual information in %f seconds\n',toc);
    end
end

%%

dI = Ihat5-Itn;
mdI = squeeze(mean(dI));
sdI = squeeze(std(dI));

figure;
errorbar(8*ks,mdI,sdI/sqrt(10));
set(gca,'XScale','log','XTick',8*ks);
xlim([0 1e4]);
xlabel('N');
ylabel('Absolute Error (nats)');

%%

pI = 100*(Ihat5-Itn)/Itn;
mpI = squeeze(mean(pI));
spI = squeeze(std(pI));

figure;
errorbar(8*ks,mpI,spI/sqrt(10));
set(gca,'XScale','log','XTick',8*ks);
xlim([0 1e4]);
xlabel('N');
ylabel('Relative Error (%)');

%%

Ihat6 = zeros(10000,1);
s = kron((1:8)',ones(125,1));
X = zeros(1000,1);
Y = zeros(1000,3);

for ii = 1:10000
    fprintf('Mixed-pair RV, iteration #%d\n',ii);
    
    tic;
    for jj = 1:8
        X(125*(jj-1)+1:125*jj) = binornd(3,jj/9,125,1);

        for kk = 1:3
            idx = find(X(125*(jj-1)+1:125*jj) == kk);
            mu = jj*ones(1,kk);
            sigma = eye(kk);
            Y(125*(jj-1)+idx,1:kk) = mvnrnd(mu,sigma,numel(idx));
        end
    end

    S = cell(1,4);
    R = cell(1,4);

    for jj = 1:4
        idx = X == jj-1;
        S{jj} = s(idx);
        R{jj} = Y(idx,1:max(1,jj-1));
    end
    fprintf('Generated random sample in %f seconds\n',toc);
    
    tic;
    Ihat6(ii) = victorMutualInformation(S,R,[1 0 0 0],true);
    fprintf('Calculated mutual information in %f seconds\n',toc);
end

%%

dI3 = zeros(10,6);

for ii = 1:6
    n = ns(ii);
    idx = randperm(10000);
    
    for jj = 1:10
        sample = Ihat6(idx(n*(jj-1)+1:n*jj),:);
        sampleMean = squeeze(mean(sample))';
        
        dI3(jj,ii) = abs(sampleMean-Itn);
    end
end

pI3 = 100*dI3/Itn;

%%

mdI = squeeze(mean(dI3));
sdI = squeeze(std(dI3));
figure;
errorbar(ns,mdI,sdI);
set(gca,'XScale','log','XTick',ns);
xlim([0 1e3]);
xlabel('N');
ylabel('Absolute Error (nats)');

%%

mpI = squeeze(mean(pI3));
spI = squeeze(std(pI3));

figure;
errorbar(ns,mpI,spI);
set(gca,'XScale','log','XTick',ns);
xlim([0 1e3]);
xlabel('N');
ylabel('Relative Error (%)');