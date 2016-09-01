wildTypeFile = 'Stim 6 4V';
fullFieldsFile = 'Stim 3 0 18BGA 8V';

%%

extractMCDSpikes(fullFieldsFile);
extractVSyncTimings(fullFieldsFile,'method','threshold','threshold',1);

responsiveCells = analyseULEDSquareResponses(fullFieldsFile,fullFieldsFile,'stim_file_20140508','blankttls','yes','figs','no','overwrite','yes','functions',{'squaresnew'},'suffixes',1);
responsiveChannels = num2str(responsiveCells(:,1));

mungeWaveClusToNex(fullFieldsFile,responsiveChannels,'responsive');

extractMCDSpikes(wildTypeFile);
extractVSyncTimings(wildTypeFile,'method','threshold','threshold',1);
mungeWaveClusToNex(wildTypeFile,responsiveChannels,'responsive');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                         %
% Simultaneously spike sort the two files %
%                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

mungeNexWaveformsToWaveClus(wildTypeFile);
mungeNexWaveformsToWaveClus(fullFieldsFile);

responsiveCells = analyseULEDSquareResponses(fullFieldsFile,fullFieldsFile,'stim_file_20140508','blankttls','yes','figs','no','overwrite','yes','functions',{'squaresnew'},'suffixes',1,'responsivecells',[],'forceclustered',true);

%%
analyseULEDSquareResponses(wildTypeFile,wildTypeFile,'','blankttls','yes','figs','no','overwrite','yes','functions',{'squaresta' 'barsta'},'suffixes',1:2,'prefixes',{'stim_file_20140710' 'stim_file_20140508'},'responsivecells',responsiveCells);