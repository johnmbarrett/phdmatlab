lambdaS = [2:5 ones(1,4); ones(1,4) 2:5];
% lambda = mean(lambdaS);

T = 1;
deltaT = 1e-3;

% Hr = lambda*T*(1-log(lambda));
Hr_s = sum(lambdaS.*T.*(1-log(lambdaS))/8,2);

% Irs = Hr - Hr_s;
% Irs2 = Irs/log(2);

%%

Hnr = zeros(2,200);
Hlr = zeros(2,200);
Dkl = zeros([size(lambdaS) 200]);

for ii = 1:200
    pii_S = exp(-lambdaS*T).*(lambdaS*T).^ii/factorial(ii);
    pii = mean(pii_S,2);
    
    for jj = 1:size(lambdaS,2)
        pii_s = pii_S(:,jj);
        Dkl(:,jj,ii) = plogq(pii_s,pii_s./pii);
    end
    
    Hnr(:,ii) = -plogp(pii);
    Hlr(:,ii) = -pii.*(sum(log(1:ii)) - ii*log(T));
end

%%

Hnr = sum(Hnr,2);
Hlr = sum(Hlr,2);
Hr = Hnr + Hlr;

Irs = (Hr - Hr_s)/log(2);

%%

HnR = zeros(200,200);
HlR = zeros(200,200);

for ii = 1:200
    for jj = 1:200
        pij = mean(exp(-sum(lambdaS,1)*T).*(lambdaS(1,:)*T).^ii.*(lambdaS(2,:)*T).^jj/(factorial(ii).*factorial(jj)),2);
        
        HnR(ii,jj) = -plogp(pij);
        HlR(ii,jj) = -pij.*(sum(log(1:ii)) + sum(log(1:jj)) - ii*log(T) - jj*log(T));
    end
end

HlR(isnan(HlR)) = 0;

%%

HnR = sum(sum(HnR));
HlR = sum(sum(HlR));
HR = HnR + HlR;

HR_s = sum(Hr_s);

IRs = (HR - HR_s)/log(2);

%%

Dkl = sum(Dkl,3);
red = mean(min(Dkl,[],1),2);

%%

unq = Irs-red;
syn = IRs-sum(unq)-red;

%%

figure;
bar([red;unq;syn;IRs])

%%

assert(false);

%%

ns = 8*ceil(logspace(0.5,4,8)/8);
ks = ns/8;

%%

Ihat = zeros(10,numel(ns),2);

for ii = 1:numel(ns)
    n = ns(ii);
    k = ks(ii);
    
    for jj = 1:10
        tic;
        s = kron((1:8)',ones(k,1));
        R = simulatePoissonProcess(lambdaS(s)',T,deltaT);
        
        Ihat(jj,ii,1) = binlessInfo(s,R,'stratification_strategy',1,'max_embed_dim',2,'min_embed_dim',1);
        Ihat(jj,ii,2) = binlessInfo(s,R,'stratification_strategy',2,'max_embed_dim',2,'min_embed_dim',1);
        toc;
    end
end