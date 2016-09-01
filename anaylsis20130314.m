%%
% recordings = initRecordings;
% locateReceptiveField(recordings(5),'chlabel',[12 13 14 15 16 17 22 23 26 27 28 31 33 34 35 41 43 44 45 46 51 53 56 58 62 64 68  71 73 75 76 77 78 82 83 84]);
batchFile(7,'cpsrh','yes','ppsrh','yes','forceclustered','yes','chlabel',[12 13 14 15 16 17 22 23 26 27 28 31 33 34 35 41 43 44 45 46 51 53 56 58 62 64 68  71 73 75 76 77 78 82 83 84]);
%%
cd ../JBWT0021;
%%
recordings = initRecordings;
recordings(6).rasterFn = getMovingBarRasterFn(3);
recordings(6).factors = {{'startX' 'startY'} 'speed' 'length' 'width'};
save('initRecordings.mat','recordings');
%%
% locateReceptiveField(recordings(3),'chlabel',[17 27 28 54 62 71]);
%%
batchFile(6,'cpsrh','yes','ppsrh','yes','forceclustered','yes','chlabel',[17 27 28 54 62 71]);
% batchFile(6,'cpsrh','yes','ppsrh','yes','forceclustered','yes','chlabel',54);