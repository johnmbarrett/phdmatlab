pws = [5 10 25 50 75 100];
nPWs = numel(pws);
reps = 10;
outfileprefix = 'stim_file_20140710';
totalTime = 0;
ISI = 60;
% fout = 1;

%% RECEPTIVE FIELDS
sizes = [1 2 4];
xs = 0:15;
ys = 0:15;

conditions = zeros(0,3);

for ii = 1:numel(sizes)
    n = 16-sizes(ii)+1;
    
    conditions = [conditions; ...
            repmat(xs(1:end-sizes(ii)+1)',n,1) ...
            kron(ys(1:end-sizes(ii)+1)',ones(n,1)) ...
            repmat(sizes(ii),n*n,1) ...
        ]; %#ok<AGROW>
end

nConditions = size(conditions,1);
conditionOrder = [];

for ii = 1:reps
    conditionOrder = [conditionOrder; conditions(randperm(nConditions),:)]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,1),'w');
pw = pws(end);

pixels = zeros(16,16,size(conditionOrder,1));

for hh = 1:size(conditionOrder,1)
    x = conditionOrder(hh,1);
    y = conditionOrder(hh,2);
    s = conditionOrder(hh,3);
    
    pixels(x+(1:s),y+(1:s),hh) = 1;
    fprintf(fout,'r %d %d %d %d %d\r\n',x,y,s,s,pw);
end

fclose(fout);

stimTime = pw*size(conditionOrder,1);
disp(stimTime/1000);
disp(stimTime/(1000*60));