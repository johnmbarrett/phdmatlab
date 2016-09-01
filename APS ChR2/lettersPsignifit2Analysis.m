letterExperiments = [56:59 61:63];
nLetterExperiments = numel(letterExperiments);

fitParams = zeros(nLetterExperiments,2,2);

% fin = fopen('psignifit2_params.txt');
% psiParams = fread(fin,[1 Inf],'char=>char');

%%

colours = 'bgrm';
za = 2*log(1/0.1-1);

for ii = 1:nLetterExperiments;
    A = importdata(sprintf('JBOG%04d_letters_psignifit_data.txt',letterExperiments(ii)));
    
    figure;
    ax = axes; %#ok<LAXES>
    hold on;
    xlim([0 32]);
    ylim([0 1]);
    
    for jj = 1:2
        scatter(ax,A(:,1),A(:,1+jj),colours(jj));
        figure;
        s = pfit(A(:,setdiff(1:4,4-jj)),'plot','shape','logistic','n_intervals',10,'fix_lambda',0,'alpha_prior','-gaussian 16, 8','beta_limits',[0 360/za]);
        fitParams(ii,:,jj) = s.params.est(1:2);
        fplot(ax,@(x) 0.1+0.9./(1+exp(-(x-s.params.est(1))/s.params.est(2))),[0 32],'Color',colours(jj+2));
    end
end

%%

midpoints = squeeze(fitParams(:,1,:));
widths = 2*log(1/0.1-1)*squeeze(fitParams(:,2,:));
lambdas = zeros(nLetterExperiments,2);

%%
% copied and pasted from aggregateAPSChR2Results.m

figure;
hold on;

colours = 'bg';
sigmoidFun = @(m,w,l,a) @(x) 0.1+(0.9-l)./(1+exp(-2.*log(1/a-1).*(x-m)./w));
inverseFun = @(m,w,l,a) @(y) m-w.*log((0.9-l)./(y-0.1)-1)./(2*log(1./a-1));

maxLogMar = ceil(10*log(160*60/5)/log(10))/10;

thresh20s = zeros(nLetterExperiments,2);
logMARs = maxLogMar*ones(nLetterExperiments,2)+0.1;

for ii = 1:2
    for jj = 1:nLetterExperiments
        sigmoid = sigmoidFun(midpoints(jj,ii),widths(jj,ii),lambdas(jj,ii),0.1);
        fplot(sigmoid,[0 20],colours(ii));
        inverse = inverseFun(midpoints(jj,ii),widths(jj,ii),lambdas(jj,ii),0.1);
        thresh20s(jj,ii) = ceil(10*log(inverse(0.2)*60)/log(10))/10;
        
        logMarRange = thresh20s(jj,ii):0.1:maxLogMar;
        
        if isempty(logMarRange)
            continue;
        end
        
        angleRange = (10.^logMarRange)/60;
        logMARs(jj,ii) = maxLogMar+0.1-0.02*sum(5*sigmoid(angleRange));
    end
end

ylim([0 1])

snellenDenominators = 20*10.^logMARs;

%%

disp(midpoints);
disp(widths);
disp(logMARs)
disp(snellenDenominators)
disp(100*mean(exp(-diff(logMARs,[],2))-1))

%%

save('letters_psignifit2_results.mat','midpoints','widths','logMARs','snellenDenominators','sigmoidFun','inverseFun');

%%

for ii = 1:nLetterExperiments
    fprintf('\t%d\t',ii);
    
    for jj = 1:2
        fprintf('\t& %4.2f\t& %4.2f\t& %3.2f\t& 20/%05u',midpoints(ii,jj),widths(ii,jj),logMARs(ii,jj),uint32(snellenDenominators(ii,jj)));
    end
    
    fprintf(' \\\\\n');
end

fprintf('\t\\hline\n\tMedian');

for ii = 1:2
    fprintf('\t& %4.2f\t& %4.2f\t& %3.2f\t& 20/%05u',median(midpoints(:,ii)),median(widths(:,ii)),median(logMARs(:,ii)),uint32(median(snellenDenominators(:,ii))));
end

fprintf(' \\\\\n');