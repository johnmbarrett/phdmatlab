%% FULL FIELD FLASHES

pws = [5 10 25 50 75 100];
nPWs = numel(pws);
reps = 20;
IPI = 2000;
outfileprefix = 'stim_file_20140212';

pwOrder = [];
for ii = 1:reps
    pwOrder = [pwOrder pws(randperm(nPWs))]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,1),'w');

for ii = 1:numel(pwOrder)
    pw = pwOrder(ii);
    
    fprintf(fout,'r 0 0 16 16 %d\r\n',pw);
    fprintf(fout,'b %d\r\n',IPI-pw);
end

fclose(fout);

%% MOVING BARS

ISI = 4000;
periods = [50 100];
nPeriods = numel(periods);

width = 2;

directions = {'we' 'nwse' 'ns' 'nesw' 'ew' 'senw' 'sn' 'swne'};
nDirections = numel(directions);
% disp(directions);
% disp(nDirections);

conditions = [repmat((1:nPeriods)',nDirections,1) kron((1:nDirections)',ones(nPeriods,1))];

reps = 10;
conditionOrder = zeros(0,2);

for ii = 1:reps
    conditionOrder = [conditionOrder; conditions(randperm(size(conditions,1)),:)]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,2),'w');

for ii = 1:size(conditionOrder,1)
    index = conditionOrder(ii,:);
    period = periods(index(1));
    direction = directions{index(2)};
%     disp(index);
%     disp(direction);
    
    fprintf(fout,'m %s %d %d\r\n',direction,period,width);
    fprintf(fout,'b %d\r\n',ISI-(16-width)*period);
end

fclose(fout);

%% FLASHING CIRCLES/RINGS

ISI = 2000;
period = 100;

centres = [4 4; 4 11; 11 4; 11 11; 8 8];
maxRadii = [2 4 8];
widths = [2 4 8];

% conditions = [repmat((1:5)',4,1) repmat(kron([1; 2],ones(5,1)),2,1) kron([1; 2],ones(10,1)); 5*ones(3,1) 3*ones(3,1) (1:3)'];
conditions = [...
    (1:5)' ones(5,1) ones(5,1); ...
    repmat((1:5)',2,1) 2*ones(10,1) kron([1;2],ones(5,1)); ...
    5*ones(3,1) 3*ones(3,1) (1:3)' ...
    ];

reps = 10;
conditionOrder = zeros(0,3);

for ii = 1:reps
    conditionOrder = [conditionOrder; conditions(randperm(size(conditions,1)),:)]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,3),'w');

for ii = 1:size(conditionOrder,1)
    index = conditionOrder(ii,:);
    centre = centres(index(1),:);
    x = centre(1);
    y = centre(2);
    radius = maxRadii(index(2));
    width = widths(index(3));
    
    fprintf(fout,'c %d %d %d %d %d %d\r\n',x,y,radius,radius,width,period);
    fprintf(fout,'b %d\r\n',ISI-period);
end

fclose(fout);