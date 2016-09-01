function board = checkBoard(width,height,x,y)
    [X,Y] = meshgrid(0:width-1,0:height-1);
    board = 255*(mod(floor(Y/y),2) == mod(floor(X/x),2));
end