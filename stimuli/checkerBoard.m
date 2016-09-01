function board = checkBoard(X,Y,x,y)
    board = meshgrid(1:X,1:Y);
    board = mod(floor(board/x),2) & mod(floor(board/y),2);
end