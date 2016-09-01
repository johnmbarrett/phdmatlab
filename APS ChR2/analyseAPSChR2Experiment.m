isFullFields = false;
isLetters = true;
firstGratings = 'luminance'; % 'contrast';

assert(isLetters);

%%

detectChR2ResponsiveCells(1+isFullFields,0,49,6,5,['./' firstGratings ' sequence.mat'],[firstGratings '_ctrl'],true,false,true,50,2);
detectChR2ResponsiveCells(2+isFullFields,0,49,6,5,'./frequency sequence.mat','frequency_ctrl',true,false,true,50,2);

detectChR2ResponsiveCells(3+isFullFields+isLetters,0,49,6,5,['./' firstGratings ' sequence.mat'],[firstGratings '_drug'],true,false,true,50,2);
detectChR2ResponsiveCells(4+isFullFields+isLetters,0,49,6,5,'./frequency sequence.mat','frequency_drug',true,false,true,50,2);

% assert(false);

%%

nCells = [487 231];
cells = arrayfun(@(n) 1:n,nCells,'UniformOutput',false);

% for ii = 1:2
    for jj = 1:1200
%         for kk = 1:nCells(ii)
%             tic;
            grating2AFCAnalysis(firstGratings,1+isFullFields,3+isFullFields+isLetters,false,25,'discrete',struct('section','main','stimulusIter',jj,'conditionIter',1:2,'cellIter',{cells}));
%             toc;
%         end
    end
% end

%%

% grating2AFCAnalysis(firstGratings,1+isFullFields,3+isFullFields+isLetters,false,25,'jacobskde',struct('section','init'));
grating2AFCAnalysis('frequency',2+isFullFields,4+isFullFields+isLetters,false,25,'jacobskde',struct('section','init'));

%%

% grating2AFCAnalysis(firstGratings,1+isFullFields,3+isFullFields+isLetters,false,25,'jacobskde',struct('section','main','stimulusIter',1:1200,'conditionIter',1:2,'cellIter',{{1:487 1:231}}));
grating2AFCAnalysis('frequency',2+isFullFields,4+isFullFields+isLetters,false,25,'jacobskde',struct('section','main','stimulusIter',1:1200,'conditionIter',1:2,'cellIter',{{1:307 1:125}}));

%%

% grating2AFCAnalysis(firstGratings,1+isFullFields,3+isFullFields+isLetters,false,25,'jacobskde',struct('section','finish'));
grating2AFCAnalysis('frequency',2+isFullFields,4+isFullFields+isLetters,false,25,'jacobskde',struct('section','finish'));

%%

grating2AFCAnalysis('frequency',2+isFullFields,4+isFullFields+isLetters,false);

%%

% contrastFrequencyGratingAnalysis('contrast',1+isFullFields,4+isFullFields,false,49);
% contrastFrequencyGratingAnalysis('frequency',2+isFullFields,5+isFullFields,false);
% 
% %%
% 
% detectChR2ResponsiveCells(4,0,17,150,5,'V:\retina\John B\phd backup\aps stimuli\frontiers paper\moving bars\sequence.mat','bars_ctrl',true,false,true);
% detectChR2ResponsiveCells(7,0,17,150,5,'V:\retina\John B\phd backup\aps stimuli\frontiers paper\moving bars\sequence.mat','bars_drug',true,false,true);
% 
% movingBarAnalysis;

%%

assert(isLetters);
detectChR2ResponsiveCells(3+isFullFields,0,5,4,5,'letters sequence.mat','letters_ctrl',true,false,true);
detectChR2ResponsiveCells(6+isFullFields,0,5,4,5,'letters sequence.mat','letters_drug',true,false,true);

%%

sloanLettersAnalysis(3+isFullFields,6+isFullFields);