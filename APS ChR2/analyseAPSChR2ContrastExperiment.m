% preStims = 0;

%%

if ~exist('skipDetection','var') || ~skipDetection
    detectChR2ResponsiveCells(1+preStims,0,9,4,5,'./lumstep sequence','lumstep_ctrl',true,false,true,50,1);
    detectChR2ResponsiveCells(2+preStims,0,26,12,5,'./constep sequence','constep_ctrl',true,false,true,25,1);
    detectChR2ResponsiveCells(3+preStims,0,7,12,5,'./revgrat sequence','revgrat_ctrl',true,false,true,50,1);

    detectChR2ResponsiveCells(4+preStims,0,9,4,5,'./lumstep sequence','lumstep_drug',true,false,true,50,1);
    detectChR2ResponsiveCells(5+preStims,0,26,12,5,'./constep sequence','constep_drug',true,false,true,25,1);
    detectChR2ResponsiveCells(6+preStims,0,7,12,5,'./revgrat sequence','revgrat_drug',true,false,true,50,1);
end

%%

luminanceStepsAnalysis(1+preStims,4+preStims,'ctrl','drug');
% contrastStepsAnalysis(2+preStims,5+preStims,'ctrl','drug');
% reversingGratingsAnalysis(3+preStims,6+preStims,'ctrl','drug');

%%

disp('This should always happen');

%%

return

%%

disp('This should never happen');

%% these blocks are only for special experiments, which is why they're after the return

detectChR2ResponsiveCells(7,0,9,4,5,'lumstep sequence','lumstep_blok',true,false,true,50,1);

%%

luminanceStepsAnalysis(4,7,'drug','blok','drug_vs_blok');

%%

detectChR2ResponsiveCells(1,0,5,4,5,'./letters sequence.mat','letters_lowc',true,false,true);
detectChR2ResponsiveCells(9,0,5,4,5,'./letters sequence.mat','letters_lowd',true,false,true);

%%

sloanLettersAnalysis(1,9,'./letters sequence.mat',[],'lowc','lowd','lowc_vs_lowd');

%%

detectChR2ResponsiveCells(2,0,9,4,5,'./lumstep sequence','lumstep_lowc',true,false,true,50,1);
detectChR2ResponsiveCells(10,0,9,4,5,'./lumstep sequence','lumstep_lowd',true,false,true,50,1);

%%

luminanceStepsAnalysis(2,10,'lowc','lowd','lowc_vs_lowd');