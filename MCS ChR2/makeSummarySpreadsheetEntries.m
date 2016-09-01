cd('D:\John B\Electrophysiology\Optogenetics\');
[~,~,recinfo] = xlsread('responsive cell recordings.xlsx','A1:H22');
recinfo(5:end+1,:) = recinfo(4:end,:);
recinfo(4,:) = {31 'SKF' 'A' 8 0 3 12 recinfo{5,end}};
recinfo{5,end} = NaN;
recinfo(8:end+1,:) = recinfo(7:end,:);
recinfo(7,:) = {33 'SKF' 'B' 8 16 19 28 NaN};

%%

experiments = vertcat(recinfo{2:end,1});
uniqueExperiments = unique(experiments);
nExperiments = numel(uniqueExperiments);
rows = 32*nExperiments;

recsummary = cell(rows,18);

%%

rowsSoFar = 0;
for ii = 2:size(recinfo)
    if strcmp(recinfo{ii,3},'A') || sum(experiments == recinfo{ii,1}) == 1
        recsummary{rowsSoFar+1,1} = sprintf('JBOG%04d',recinfo{ii,1});
        recsummary{rowsSoFar+1,2} = 'Frontiers paper';
        recsummary{rowsSoFar+1,4} = 'ChR2rd1';
        recsummary{rowsSoFar+1,7} = recinfo{ii,8};
        isSecondRetina = false;
    else
        isSecondRetina = true;
    end
    
    recsummary{rowsSoFar+1,8} = recinfo{ii,3};
    recsummary{rowsSoFar+1,11} = 'std';
    recsummary{rowsSoFar+1,12} = 'Inv';
    
    isSKF = strcmp(recinfo{ii,2},'SKF');
    
    middleVolt = recinfo{ii,4};
    offset = recinfo{ii,5};
    
    recsummary{rowsSoFar+1,13} = sprintf('Spont Control %d',1+isSecondRetina);
    recsummary{rowsSoFar+1,16} = 'N/A';
    recsummary{rowsSoFar+1,17} = 'Spontaneous Activity';
    
    for jj = 1:5
        recsummary{rowsSoFar+jj+1,13} = sprintf('Stim %d 0 %s %dV',jj+offset,recinfo{ii,2},5+jj);
        recsummary{rowsSoFar+jj+1,16} = 'stim_file_20140508_1.bat';
        recsummary{rowsSoFar+jj+1,17} = 'Full-field flashes';
        recsummary{rowsSoFar+jj+1,18} = 5+jj-4*(jj==6);
    end
    
    recsummary{rowsSoFar+7,13} = sprintf('Stim %d 0 %s %dV',6+offset,recinfo{ii,2},middleVolt);
    recsummary{rowsSoFar+7,16} = 'stim_file_20140710_2.bat';
    recsummary{rowsSoFar+7,17} = 'Flashing squares, moving bars';
    recsummary{rowsSoFar+7,18} = middleVolt;
    
    for jj = 1:3+isSKF
        recsummary{rowsSoFar+2*(jj-1)+8,11} = sprintf('std + %duM %s',10*2^(jj-1),recinfo{ii,2});
        recsummary{rowsSoFar+2*(jj-1)+8,13} = sprintf('Spont %duM %s',10*2^(jj-1),recinfo{ii,2});
        recsummary{rowsSoFar+2*(jj-1)+8,16} = 'N/A';
        recsummary{rowsSoFar+2*(jj-1)+8,17} = 'Spontaneous Activity';
        
        recsummary{rowsSoFar+2*(jj-1)+9,13} = sprintf('Stim %d 0 %s %dV',jj+6+offset,recinfo{ii,2},middleVolt);
        recsummary{rowsSoFar+2*(jj-1)+9,16} = 'stim_file_20140508_1.bat';
        recsummary{rowsSoFar+2*(jj-1)+9,17} = 'Full-field flashes';
        recsummary{rowsSoFar+2*(jj-1)+9,18} = middleVolt;
    end
    
    recsummary{rowsSoFar+2*isSKF+14,11} = sprintf('std + %duM %s',10*2^(3+isSKF),recinfo{ii,2});
    recsummary{rowsSoFar+2*isSKF+14,13} = sprintf('Spont %duM %s',10*2^(3+isSKF),recinfo{ii,2});
    recsummary{rowsSoFar+2*isSKF+14,16} = 'N/A';
    recsummary{rowsSoFar+2*isSKF+14,17} = 'Spontaneous Activity';
    
    for jj = 1:5
        recsummary{rowsSoFar+jj+2*isSKF+14,13} = sprintf('Stim %d %d %s %dV',jj+9+isSKF+offset,80+80*isSKF,recinfo{ii,2},5+jj);
        recsummary{rowsSoFar+jj+2*isSKF+14,16} = 'stim_file_20140508_1.bat';
        recsummary{rowsSoFar+jj+2*isSKF+14,17} = 'Full-field flashes';
        recsummary{rowsSoFar+jj+2*isSKF+14,18} = 5+jj-4*(jj==6);
    end
    
    recsummary{rowsSoFar+20+2*isSKF,13} = sprintf('Stim %d %d %s %dV',15+isSKF+offset,80+80*isSKF,recinfo{ii,2},middleVolt);
    recsummary{rowsSoFar+20+2*isSKF,16} = 'stim_file_20140710_2.bat';
    recsummary{rowsSoFar+20+2*isSKF,17} = 'Flashing squares, moving bars';
    recsummary{rowsSoFar+20+2*isSKF,18} = middleVolt;
    
    recsummary{rowsSoFar+21+2*isSKF,13} = sprintf('Spont Washout %d',1+isSecondRetina);
    recsummary{rowsSoFar+21+2*isSKF,16} = 'N/A';
    recsummary{rowsSoFar+21+2*isSKF,17} = 'Spontaneous Activity';
    
    recsummary{rowsSoFar+22+2*isSKF,13} = sprintf('Stim %d 0 %s %dV',16+offset,recinfo{ii,2},middleVolt);
    recsummary{rowsSoFar+22+2*isSKF,16} = 'stim_file_20140508_1.bat';
    recsummary{rowsSoFar+22+2*isSKF,17} = 'Full-field flashes';
    recsummary{rowsSoFar+22+2*isSKF,18} = middleVolt;
    
    rowsSoFar = rowsSoFar + 6 + 16 + 2*isSKF;
end

%%

xlswrite('frontiers_paper_experiments_summary.xlsx',recsummary);