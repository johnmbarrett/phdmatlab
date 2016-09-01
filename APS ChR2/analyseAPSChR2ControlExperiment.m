detectChR2ResponsiveCells(1,0,49,6,5,'./frequency sequence.mat','frequency_ctrl',true,false,true,50,2);
detectChR2ResponsiveCells(4,0,49,6,5,'./frequency sequence.mat','frequency_drug',true,false,true,50,2);
detectChR2ResponsiveCells(9,0,49,6,5,'./frequency sequence.mat','frequency_high',true,false,true,50,2);

detectChR2ResponsiveCells(2,0,9,4,5,'./lumstep sequence','lumstep_ctrl',true,false,true,50,1);
detectChR2ResponsiveCells(5,0,9,4,5,'./lumstep sequence','lumstep_drug',true,false,true,50,1);
detectChR2ResponsiveCells(8,0,9,4,5,'./lumstep sequence','lumstep_high',true,false,true,50,1);

detectChR2ResponsiveCells(3,0,5,4,5,'./letters sequence.mat','letters_ctrl',true,false,true);
detectChR2ResponsiveCells(6,0,5,4,5,'./letters sequence.mat','letters_drug',true,false,true);
detectChR2ResponsiveCells(7,0,5,4,5,'./letters sequence.mat','letters_high',true,false,true);

%%

grating2AFCAnalysis('frequency',1,4,false,25,'discrete',struct('section','none'),'ctrl','drug','ctrl_vs_drug');
grating2AFCAnalysis('frequency',4,9,false,25,'discrete',struct('section','none'),'drug','high','drug_vs_high');

%%

sloanLettersAnalysis(3,6,'./letters sequence.mat',[],'ctrl','drug','ctrl_vs_drug');
sloanLettersAnalysis(6,7,'./letters sequence.mat',[],'drug','high','drug_vs_high');

%%

luminanceStepsAnalysis(2,5,'ctrl','drug','ctrl_vs_drug');
luminanceStepsAnalysis(5,8,'drug','high','drug_vs_high');