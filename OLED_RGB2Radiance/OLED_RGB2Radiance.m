%% p_RGB2ImageRadianceMap
%    
% Tutorial for how to generate an image radiance map for any given bitmap
%
% This tutorial illustrates how to use the isetBio functions to calculate
% the radiance image emitted when a calibrated display renders a bitmap.
%
% The display calibration data are stored in a Matlab file that includes 
% 1. the spectral power distribution of display color primaries
% 2. a lut mapping display rgb values to luminance intensity, or alternatively
%    the exponent of a display gamma function 
%
% Examples of display calibration data can be found in isetBio/isettools/data/displays
%
%
% isetBio Team 2013 (H.J, B.W and J.F)

%% Initialize clears the workspace and starts a new version of isetBio
s_initISET;

%% Create scene from file
%  Select the rgb image and the display data files
%  We have several rgb files (tif and jpg) stored in isetBio/isettools/data/images/rgb 
imgFileName = 'macbeth.tif'; 

% We have other display calibration data stored in isetBio/isettools/data/displays
% Here are the calibration files for the Sony BVM-F250 and PVM-2541 from
% the measurements described in Cooper, Jiang, Vildavski, Farrell, &
% Norcia, JOV, 2013
dispBVMFile = 'OLED-SonyBVM-F250.mat';
dispPVMFile = 'OLED-SonyPVM-2541.mat';

%  Check existence
if ~exist(imgFileName,'file'), error('Image file not found'); end
if ~exist(dispBVMFile,'file'), error('BVM Display file not found.'); end
if ~exist(dispPVMFile,'file'), error('PVM Display file not found.'); end

%%  Create scene from file
% Scenes are stored in isetBio as structures that include
%       a matrix of scene data represented as photons with 32 bits of precision by default
%       a vector with the wavelengths that are used in the spectral representation (400:10:700 by default, but can be changed) 
%       a vector with illuminant data (this allows users to calculate reflectance)
%               in the case where the scene is an image displayed on a monitor, we set
%               the illuminant to be the white point of the monitor.
% It is possible to set and get properties of the scene using the functions sceneGet and sceneSet
% Use "Plot" and "Analyze" pull down menus to view summaries of the
% spectrum and luminance

%%  Scene on the Sony OLED BVM
sceneB = sceneFromFile(imgFileName,'rgb',[],dispBVMFile);
sceneB = sceneSet(sceneB,'name','Scene on BVM');
vcAddAndSelectObject(sceneB); sceneWindow;

%%  Scene on the Sony OLED PVM
sceneP = sceneFromFile(imgFileName,'rgb',[],dispPVMFile);
sceneP = sceneSet(sceneP,'name','Scene on PVM');
vcAddAndSelectObject('scene',sceneP); sceneWindow;

%%  Scene on CRT
sceneC = sceneFromFile(imgFileName,'rgb',[],'CRT-Dell');
sceneC = sceneSet(sceneC,'name','Scene on CRT-Dell');
vcAddAndSelectObject('scene',sceneC); sceneWindow;

%%  Scene on other CRT
sceneC = sceneFromFile(imgFileName,'rgb',[],'LCD-Dell');
sceneC = sceneSet(sceneC,'name','Scene on LCD-Dell');
vcAddAndSelectObject('scene',sceneC); sceneWindow;

%% Compare the four images
imageMultiview('scene', 1:4, true);

%% Compare the gamuts of the four monitors
d = displayCreate('CRT-Dell');
displayPlot(d,'gamut');
title('Dell CRT')

d = displayCreate('LCD-Dell');
displayPlot(d,'gamut');
title('Dell LCD')

d = displayCreate('OLED-SonyBVM-F250');
displayPlot(d,'gamut');
title('Sony OLED (BVM)')

d = displayCreate('OLED-SonyPVM-2541');
displayPlot(d,'gamut');
title('Sony OLED (PVM)')


%% End