function cross = calibrationCross(height,width,thickness)
    cross = 255*ones(height,width);
    
    midHeight = floor((height+1)/2);
    midWidth = floor((width+1)/2);
    halfThickness = ceil((thickness-1)/2);
    
    indices = -halfThickness+(1-mod(thickness,2)):halfThickness;
    blackCols = indices+midWidth;
    blackRows = indices+midHeight; 
    cross(:,blackCols) = zeros(height,length(blackCols));
    cross(blackRows,:) = zeros(length(blackRows),width);
end