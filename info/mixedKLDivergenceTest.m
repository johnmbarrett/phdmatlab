lambda = 1:8;
p = lambda/9;

mu = 4.5;
q = 0.5;

K = p.*(log(p)-log(q)) + (1-p).*(log(1-p)-log(1-q)) + ...
    p.*(log(lambda)-log(mu)+mu./lambda-1);

ns = ceil(logspace(0.5,3.5,7));

%%

Khat = zeros(10,7,8);

for jj = 1:7
    n = ns(jj);
    
    for hh = 1:10
        for ii = 1:8
            tic;
            xp = binornd(1,p(ii),1000,1);
            yp = Inf(size(xp));
            yp(xp == 1) = exprnd(1/lambda(ii),sum(xp == 1),1);
            
            xq = binornd(1,q,1000,1);
            yq = Inf(size(xq));
            yq(xq == 1) = exprnd(1/mu,sum(xq == 1),1);
            Q = {yq(xq == 0) yq(xq == 1)};

            Khat(hh,jj,ii) = mixedKLDivergence({yp(xp == 0) yp(xp == 1)},Q,[true false],1,'both');
        end
        toc;
    end
end

%%

K = repmat(reshape(K',[1 1 8]),10,7);
 
%%

mdK = squeeze(mean(abs(K-Khat)));
sdK = squeeze(std(abs(K-Khat)));

%%

figure;
hold on;

errorbar(repmat(ns',1,8),mdK,sdK);

set(gca,'XScale','log');

xlabel('# Samples');
ylabel('Absolute error (nats)');
 
%%

mpK = squeeze(mean(100*abs(K-Khat)./Khat));
spK = squeeze(std(100*abs(K-Khat)./Khat));

%%

figure;
hold on;

errorbar(repmat(ns',1,8),mpK,spK);

set(gca,'XScale','log');

xlabel('# Samples');
ylabel('% Error');

%%

Kbar = zeros(10,8);

for hh = 1:10
    tic;
    xq = binornd(1,q,1000,1);
    yq = Inf(size(xq));
    yq(xq == 1) = exprnd(1/mu,sum(xq == 1),1);
    Q = {yq(xq == 0) yq(xq == 1)};
    
    for ii = 1:8
        xp = binornd(1,p(ii),125,1);
        yp = Inf(size(xp));
        yp(xp == 1) = exprnd(1/lambda(ii),sum(xp == 1),1);

        Kbar(hh,ii) = mixedKLDivergence({yp(xp == 0) yp(xp == 1)},Q,[true false],1,'both');
    end
    toc;
end

%%

K = repmat(K,10,1);
 
%%

mdK = squeeze(mean(abs(K-Kbar)));
sdK = squeeze(std(abs(K-Kbar)));

%%

figure;
hold on;

errorbar(1:8,mdK,sdK);

set(gca,'XScale','log');

xlabel('# Samples');
ylabel('Absolute error (nats)');
 
%%

mpK = squeeze(mean(100*abs(K-Kbar)./Kbar));
spK = squeeze(std(100*abs(K-Kbar)./Kbar));

%%

figure;
hold on;

errorbar(1:8,mpK,spK);

set(gca,'XScale','log');

xlabel('# Samples');
ylabel('% Error');

%%

ns = ceil(logspace(0.25,3.75,15));

dK = zeros(10000,15,8);

for hh = 1:10000
    tic;
    for ii = 1:15
        idx = randperm(10000);
        dK(hh,ii,:) = K-mean(Khat(idx(1:ns(ii)),:));
    end
    toc;
end

%%

mdK = squeeze(mean(dK));
mpK = 100*mdK./repmat(K,15,1);

figure;
plot(ns,mdK);
set(gca,'XScale','log');
xlabel('# samples');
ylabel('Absolute Error (nats)');

figure;
plot(ns,mpK);
set(gca,'XScale','log');
xlabel('# samples');
ylabel('% Error');