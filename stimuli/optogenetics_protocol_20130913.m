pws = [5 10 25 50 75 100];
nPWs = numel(pws);
reps = 10;
outfileprefix = 'stim_file_20130917';
% fout = 1;

%% STIMULATE SINGLE LEDs
pwOrder = repmat(1:nPWs,1,reps);
pwOrder = pwOrder(randperm(nPWs*reps));
quadrants = [0 0; 1 1; 0 1; 1 0];
mask = 3;

fout = fopen(sprintf('%s_%d.txt',outfileprefix,1),'w');

for hh = 1:nPWs*reps
    pw = pws(pwOrder(hh));
    pd = pws(end)-pw;
    
    for ii = 1:16
        for jj = 1:16
            ledIndex = sub2ind([16 16],jj,ii);
            x = 0;
            y = 0;

            for kk = 1:4
                bits = 2*(kk-1);
                quadIndex = bitshift(bitand(ledIndex-1,bitshift(mask,bits)),-bits);

                quad = quadrants(quadIndex+1,:);
                shift = 2^(4-kk);

                x = x + shift*quad(1);
                y = y + shift*quad(2);
            end

            fprintf(fout,'r %d %d 1 1 %d\r\n',x,y,pw);
            
            if pd > 0
                fprintf(fout,'b %d\r\n',pd);
            end
        end
    end
end

fclose(fout);

%% RESPONSE VS. PW/# LEDs

blocks = ...
   {0 0 16 16; ...
    0 0 8 8; ...
    0 8 8 8; ...
    8 0 8 8; ...
    8 8 8 8; ...
    0 0 4 4; ...
    0 4 4 4; ...
    0 8 4 4; ...
    0 12 4 4; ...
    4 0 4 4; ...
    4 4 4 4; ...
    4 8 4 4; ...
    4 12 4 4; ...
    8 0 4 4; ...
    8 4 4 4; ...
    8 8 4 4; ...
    8 12 4 4; ...
    12 0 4 4; ...
    12 4 4 4; ...
    12 8 4 4; ...
    12 12 4 4};
    
nBlocks = size(blocks,1);
IPI = 2000;

conditions = [repmat((1:nPWs)',nBlocks,1) kron((1:nBlocks)',ones(nPWs,1))];
conditionOrder = zeros(0,2);

for ii = 1:reps
    conditionOrder = [conditionOrder; conditions(randperm(size(conditions,1)),:)]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,2),'w');

for ii = 1:size(conditionOrder,1)
    index = conditionOrder(ii,:);
    pw = pws(index(1));
    block = blocks(index(2),:);
    
    fprintf(fout,'r %d %d %d %d %d\r\n',block{:},pw);
    fprintf(fout,'b %d\r\n',IPI-pw);
end

fclose(fout);

%% RESPONSE VS. PW/FREQ

% pws = 10;
freqs = [1 2 5 10 20 50];
nFreqs = numel(freqs);
trainLength = 5000;
ISI = 5000;

% conditions = [repmat((1:nPWs)',nFreqs,1) kron((1:nFreqs)',ones(nPWs,1))];
conditions = [repmat((2:4)',4,1) kron((1:4)',ones(3,1)); 2 5; 3 5; 2 6];
conditionOrder = zeros(0,2);

for ii = 1:reps
    conditionOrder = [conditionOrder; conditions(randperm(size(conditions,1)),:)]; %#ok<AGROW>
end

%%

fout = fopen(sprintf('%s_%d.txt',outfileprefix,3),'w');

for ii = 1:size(conditionOrder,1)
    index = conditionOrder(ii,:);
    pw = pws(index(1));
    freq = freqs(index(2));
    pd = 1000/freq;
    np = trainLength/pd;
    
    for jj = 1:np
        fprintf(fout,'r 0 0 16 16 %d\r\n',pw);
        fprintf(fout,'b %d\r\n',pd-pw);
    end
    
    fprintf(fout,'b %d\r\n',ISI);
end

fclose(fout);

% pws = 100;

%% SINGLE PULSES VS. MULTIPLE SHORT PULSES

ISI = 2000;
conditions = {...
	5 ...
	10 ...
	[5 2 5] ...
	[5 5 5] ...
	20 ...
	[10 2 10] ...
	[10 5 10] ...
	[10 10 10] ...
	[5 2 5 2 5 2 5] ...
	[5 5 5 5 5 5 5] ...
	40 ...
	[20 5 20] ...
	[20 10 20] ...
	[15 5 15 5 15] ...
	[15 10 15 10 15] ...
	[10 2 10 2 10 2 10] ...
	[10 5 10 5 10 5 10] ...
	[10 10 10 10 10 10 10] ...
	80 ...
	[40 5 40] ...
	[40 10 40] ...
	[20 5 20 5 20 5 20] ...
	[20 10 20 10 20 10 20] ...
	[10 2 10 2 10 2 10 2 10 2 10 2 10 2 10] ...
	[10 5 10 5 10 5 10 5 10 5 10 5 10 5 10] ...
	[10 10 10 10 10 10 10 10 10 10 10 10 10 10 10]};

conditionIndices = (1:numel(conditions))';
conditionOrder = [];

for ii = 1:reps
    conditionOrder = [conditionOrder; conditionIndices(randperm(numel(conditions)))]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,4),'w');

for ii = 1:size(conditionOrder,1)
    pulseTrain = conditions{conditionOrder(ii)};
    trainDuration = sum(pulseTrain);
    
    for jj = 1:numel(pulseTrain)
        pw = pulseTrain(jj);
        
        if mod(jj,2) == 0
            fprintf(fout,'b %d\r\n',pw);
        else
            fprintf(fout,'r 0 0 16 16 %d\r\n',pw);
        end
    end
    
    fprintf(fout,'b %d\r\n',ISI-trainDuration);
end

fclose(fout);

%% WHITE NOISE

% this isn't working properly and I don't understand why

pw = 100; % for some reason gifs are sped up by a factor of 5???
nIntensities = 10; % minimum frame delay in a GIF is 1/100s, i.e. 10ms
filename = 'uled_wn.gif';
duration = 15*1000; % 15*60*1000;
nFrames = duration/pw;
nSubframes = nFrames*nIntensities;
framesPerIntensity = nFrames/nIntensities;

noise = 0:nIntensities-1; % don't use max intensity so each pixel is off for at least 10ms of every 100
noise = reshape(noise,[1 1 nIntensities]);
noise = repmat(noise,[16 16 framesPerIntensity]);

rng('shuffle');

for ii = 1:16
    for jj = 1:16
        noise(ii,jj,:) = noise(ii,jj,randperm(nFrames));
    end
end

% plot(squeeze(mean(mean(noise/nIntensities))))
% ylim([0 1]);

fout = fopen(sprintf('%s_%d.txt',outfileprefix,5),'w');

for ii = 1:nFrames
    frame = noise(:,:,ii); 
    subframes = zeros(16,16,nIntensities);
    
    for jj = 1:16
        for kk = 1:16
            maxDelay = nIntensities-frame(jj,kk)+1;
            delay = randi(maxDelay)-1;
            subframes(jj,kk,(1:frame(jj,kk)-1)+delay) = 255;
        end
    end
    
%     for jj = 1:nIntensities-1
%         indices = find(frame == jj);
%         [x,y] = ind2sub([16 16],indices(:));
%         nPixels = numel(indices);
%         nPossibleStarts = nIntensities-jj+1;
%         pixelsPerStart = nPixels/nPossibleStarts;
%         
%         if pixelsPerStart >= 1.5
%             cpps = ceil(pixelsPerStart);
%             
%             for kk = 1:nPossibleStarts
%                 subindices = cpps*(kk-1)+1:min(cpps*kk,nPixels);
%                 timeIndex = (1:jj)+kk-1;
%                 
%                 for ll = subindices
%                     subframes(x(ll),y(ll),timeIndex) = 255;
%                 end
%             end
%         elseif pixelsPerStart <= 0.5
%             period = floor(1/pixelsPerStart);
%             
%             for kk = 1:nPixels
%                 subframes(x(kk),y(kk),(1:jj)+period*(kk-1)) = 255;
%             end
%         else
%             for kk = 1:nPixels
%                 subframes(x(kk),y(kk),kk+(1:jj)-1) = 255;
%             end
%         end
%     end
%     
%     indices = find(frame == nIntensities);
%     [x,y] = ind2sub([16 16],indices(:));
%     
%     for jj = 1:numel(indices)
%         subframes(x(jj),y(jj),:) = 255;
%     end

%     plot(squeeze(mean(mean(subframes/255))));
%     ylim([0 1]);
%     surf(squeeze(mean(subframes/255,3)));
%     caxis([0 1]);
%     view(2);
    
    if ii == 1
        extraOptions = {'LoopCount' Inf};
    else
        extraOptions = {'WriteMode' 'append'};
    end
    
    imwrite(reshape(subframes,[16 16 1 nIntensities]),filename,'gif','DelayTime',(pw/nIntensities)/1000,extraOptions{:});
end

fprintf(fout,'i %d V:\\retina\\John B\\phd backup\\matlab\\stimuli\\%s\r\n',pw,filename);

fclose(fout);

%% MOVING BARS

ISI = 4000;
period = 100;

widths = 1:4;
nWidths = numel(widths);

directions = {'we' 'nwse' 'ns' 'nesw' 'ew' 'senw' 'sn' 'swne'};
nDirections = numel(directions);
% disp(directions);
% disp(nDirections);

conditions = [repmat((1:nWidths)',nDirections,1) kron((1:nDirections)',ones(nWidths,1))];

conditionOrder = zeros(0,2);

for ii = 1:reps
    conditionOrder = [conditionOrder; conditions(randperm(size(conditions,1)),:)]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,6),'w');

for ii = 1:size(conditionOrder,1)
    index = conditionOrder(ii,:);
    width = widths(index(1));
    direction = directions{index(2)};
%     disp(index);
%     disp(direction);
    
    fprintf(fout,'m %s %d %d\r\n',direction,period,width);
    fprintf(fout,'b %d\r\n',ISI-(16-width)*period);
end

fclose(fout);

%% LETTERS & NUMBERS

ISI = 2000;
stimuli = char(['0'+(0:9) 'A'+(0:25)]);
nStimuli = numel(stimuli);

conditions = [repmat((1:nPWs)',nStimuli,1) kron((1:nStimuli)',ones(nPWs,1))];
conditionOrder = zeros(0,2);

for ii = 1:reps
    conditionOrder = [conditionOrder; conditions(randperm(size(conditions,1)),:)]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,7),'w');

for ii = 1:size(conditionOrder,1)
    index = conditionOrder(ii,:);
    pw = pws(index(1));
    stimulus = stimuli(index(2));
    
    fprintf(fout,'i %d V:\\retina\\John B\\phd backup\\matlab\\stimuli\\uled letters\\%s.gif\r\n',pw,stimulus);
    fprintf(fout,'b %d\r\n',ISI-pw);
end

fclose(fout);