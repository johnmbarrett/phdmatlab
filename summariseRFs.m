% srrfs = dir('srrf*mat');
% srrfs = {srrfs.name};
% data = [];
% 
% for ii = 1:numel(srrfs)
%     srrf = srrfs{ii};
%     
%     try
%         load(srrf)
%     catch err
%         continue;
%     end
%     
%     on = pmax < 0.05;
%     off = pmin < 0.05;
%     
%     if ~on && ~off
%         continue;
%     end
%     
%     [~,~,~,~,chcl] = regexp(srrf,'channel_([0-9]+)_cluster_([0-9]+)');
%     channel = str2double(chcl{1}{1});
%     cluster = str2double(chcl{1}{2});
%     
%     if on
%         ton = 1000*(31-centreTOn)/30;
%         rfdiamon = 200*mean([sdXOn sdYOn]);
%         
%         if sdXOn < sdYOn
%             eccon = sqrt(1-(sdXOn/sdYOn)^2);
%         else
%             eccon = sqrt(1-(sdYOn/sdXOn)^2);
%         end
%     else
%         ton = nan;
%         rfdiamon = nan;
%         eccon = nan;
%     end
%     
%     if off
%         toff = 1000*(31-centreTOff)/30;
%         rfdiamoff = 200*mean([sdXOff sdYOff]);
%         
%         if sdXOff < sdYOff
%             eccoff = sqrt(1-(sdXOff/sdYOff)^2);
%         else
%             eccoff = sqrt(1-(sdYOff/sdXOff)^2);
%         end
%     else
%         toff = nan;
%         rfdiamoff = nan;
%         eccoff = nan;
%     end
%     
%     if on && off
%         onoffratio = max(max(RF(:,:,round(centreTOn))))/abs(min(min(RF(:,:,round(centreTOff)))));
%     else
%         onoffratio = NaN;
%     end
%     
%     d = [channel cluster pmax pmin ton rfdiamon eccon toff rfdiamoff eccoff onoffratio];
%     
%     if isempty(data)
%         data = d;
%     else
%         data = [data; d]; %#ok<AGROW>
%     end
% end

%%

for ii = 5:7
    x = [dark(isfinite(dark(:,ii)),ii); dark(isfinite(dark(:,ii+3)),ii+3); light(isfinite(light(:,ii)),ii); light(isfinite(light(:,ii+3)),ii+3)];
    f = [repmat({'ON'},sum(isfinite(dark(:,ii))),1); repmat({'OFF'},sum(isfinite(dark(:,ii+3))),1); repmat({'ON'},sum(isfinite(light(:,ii))),1); repmat({'OFF'},sum(isfinite(light(:,ii+3))),1)];
    g = [repmat({'Dark'},sum(isfinite(dark(:,ii)))+sum(isfinite(dark(:,ii+3))),1); repmat({'Light'},sum(isfinite(light(:,ii)))+sum(isfinite(light(:,ii+3))),1)];
    [~,~,stats] = anovan(x,{f g},'model','interaction','random',[],'varnames',{'Light Condition' 'On/Off'});
%     figure
%     hist(stats.resid);
    [~,p] = HartigansDipSignifTest(stats.resid,1000);
    disp(p);
%     figure
%     multcompare(stats,'dimension',1);
%     figure
%     multcompare(stats,'dimension',2);
%     figure
%     multcompare(stats,'dimension',[1 2]);
end

%%

[p,~,stats] = ranksum(dark(isfinite(dark(:,11)),11), light(isfinite(light(:,11)),11));
disp(p);
disp(stats.ranksum);
disp(stats.zval);

%%
figure
colours = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
ylabels = {'Peak Time Before Spike/ms', 'RF Diameter/{\mu}m', 'Eccentricity', 'ON/OFF Ratio'};

for ii = 1:3
    datas = {dark(isfinite(dark(:,ii+4)),ii+4); dark(isfinite(dark(:,ii+7)),ii+7); light(isfinite(light(:,ii+4)),ii+4); light(isfinite(light(:,ii+7)),ii+7)};
    subplot('Position',[0.5*(1-mod(ii,2))+0.085 1-floor((ii-1)/2)*0.5-0.425 0.4 0.4]); %2,2,ii);
    hold on;
    
    m = zeros(4,1);
    e = zeros(4,1);
    
    for jj = 1:4
        mu = mean(datas{jj});
        sem = std(datas{jj})/sqrt(numel(datas{jj}));
        
        m(jj) = mu;
        e(jj) = 2*sem;
        
        bar(jj,mu,'FaceColor',colours(jj,:));
    end
    
    errorbar(1:4,m,e,'Color','k','LineStyle','none');
    xlim([0.5 4.5]);
    set(gca,'XTick',1:4,'XTickLabel',{'Dark/ON','Dark/OFF','Light/ON','Light/OFF'});
    ylabel(ylabels{ii});
end

subplot('Position',[0.5*(1-mod(4,2))+0.085 1-floor((4-1)/2)*0.5-0.425 0.4 0.4]); %2,2,4);
x1 = dark(isfinite(dark(:,11)),11);
x2 = light(isfinite(light(:,11)),11);
mu = [mean(x1); mean(x2)];
sem = [std(x1)/numel(x1); std(x2)/numel(x2)];
hold on;
bar(1,mu(1),'FaceColor',colours(1,:));
bar(2,mu(2),'FaceColor',colours(2,:));
errorbar(1:2,mu,2*sem,'Color','k','LineStyle','none');
set(gca,'XTick',1:2,'XTickLabel',{'Dark','Light'});
ylabel(ylabels{4});
xlim([0.5 2.5]);