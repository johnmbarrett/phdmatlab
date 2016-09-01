pws = [5 10 25 50 75 100];
nPWs = numel(pws);
reps = 10;
outfileprefix = 'stim_file_20140508';
totalTime = 0;
ISI = 60;
% fout = 1;

%% FULL FIELD FLASHES

IPI = 2000;

pwOrder = [];
for ii = 1:2*reps
    pwOrder = [pwOrder; pws(randperm(nPWs))']; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,1),'w');

for ii = 1:numel(pwOrder)
    pw = pwOrder(ii);
    
    fprintf(fout,'r 0 0 16 16 %d\r\n',pw);
    fprintf(fout,'b %d\r\n',IPI-pw);
end

fclose(fout);

stimTime = size(pwOrder,1)*IPI/1000;
disp(stimTime/60);
totalTime = totalTime + stimTime + ISI;

%% STIMULATE SINGLE LEDs
sizes = [1 2];
conditions = [repmat(1:nPWs,1,2); kron(sizes,ones(1,nPWs))]';
conditionOrder = [];

for ii = 1:reps
    conditionOrder = [conditionOrder; conditions(randperm(size(conditions,1)),:)]; %#ok<AGROW>
end

quadrants = [0 0; 1 1; 0 1; 1 0];
mask = 3;

fout = fopen(sprintf('%s_%d.txt',outfileprefix,4),'w');

I = zeros(16,16,1,sum((16./sizes).^2)*nPWs*reps);

nn = 0;
for hh = 1:size(conditionOrder,1)
    pw = pws(conditionOrder(hh,1));
    s = conditionOrder(hh,2);
    pd = pws(end)-pw;
    
    for ii = 1:16/s
        for jj = 1:16/s
            ledIndex = sub2ind([16 16],jj,ii);
            x = 0;
            y = 0;

            for kk = 1:5-s
                bits = 2*(kk-1);
                quadIndex = bitshift(bitand(ledIndex-1,bitshift(mask,bits)),-bits);

                quad = quadrants(quadIndex+1,:);
                shift = 2^(4-kk);

                x = x + shift*quad(1);
                y = y + shift*quad(2);
            end

            nn = nn + 1;
%             I(x+(1:s),y+(1:s),1,nn) = 254;
            fprintf(fout,'r %d %d %d %d %d\r\n',x,y,s,s,pw);
            
            if pd > 0
                fprintf(fout,'b %d\r\n',pd);
            end
        end
    end
end

fclose(fout);

stimTime = nn*pws(end)/1000;
disp(stimTime/60);
totalTime = totalTime + stimTime + ISI;

%% RESPONSE VS. PW/# LEDs

blocks = [ ...
    0 0 16 16; ...
    [kron([0; 8],ones(2,1)) repmat([0; 8],2,1) 8*ones(4,2)]; ...
    [kron((0:4:12)',ones(4,1)) repmat((0:4:12)',4,1) 4*ones(16,2)]; ...
    ];
    
nBlocks = size(blocks,1);
IPI = 2000;

conditions = [repmat((1:nPWs)',nBlocks,1) kron((1:nBlocks)',ones(nPWs,1))];
conditionOrder = zeros(0,2);

for ii = 1:reps
    conditionOrder = [conditionOrder; conditions(randperm(size(conditions,1)),:)]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,5),'w');

for ii = 1:size(conditionOrder,1)
    index = conditionOrder(ii,:);
    pw = pws(index(1));
    block = num2cell(blocks(index(2),:));
    
    fprintf(fout,'r %d %d %d %d %d\r\n',block{:},pw);
    fprintf(fout,'b %d\r\n',IPI-pw);
end

fclose(fout);

nTrials = nBlocks*nPWs*reps;
stimTime = nTrials*IPI/1000;

disp(stimTime/60);
totalTime = totalTime + stimTime + ISI;

%% WHITE NOISE

filename = 'wn025.gif';

fout = fopen(sprintf('%s_%d.txt',outfileprefix,6),'w');
fprintf(fout,'g V:\\retina\\John B\\phd backup\\matlab\\stimuli\\%s\r\n',pw,filename);

fclose(fout);

% gif = imread(filename);

stimTime = size(gif,4)*0.01;

disp(stimTime/60);
totalTime = totalTime + stimTime + ISI;

%% MOVING BARS

IPI = 4000;
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
    fprintf(fout,'b %d\r\n',IPI-(16-width)*period);
end

fclose(fout);

stimTime = IPI*size(conditionOrder,1)/1000;

disp(stimTime/60);
totalTime = totalTime + stimTime + ISI;

%% FLASHING CIRCLES/RINGS

IPI = 2000;
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
    fprintf(fout,'b %d\r\n',IPI-period);
end

fclose(fout);

stimTime = IPI*size(conditionOrder,1)/1000;

disp(stimTime/60);
totalTime = totalTime + stimTime;

%%

% disp(totalTime);
disp(totalTime/60);
disp(totalTime/3600);