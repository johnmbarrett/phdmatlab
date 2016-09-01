getPixels = cell(4,1);

angles = 0:pi/2:3*pi/2;

for ii = 1:4
    getPixels{ii} = squareWave(255,64,32,angles(ii),15*60);
end

stimulate(getPixels,[0 0 640 480],[1920 0 1920+640 480],[],NaN,0,60,NaN,false,false,false,cell(size(getPixels)));