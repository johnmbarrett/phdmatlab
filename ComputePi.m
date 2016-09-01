function [ output ] = ComputePi(count)
    % Check we have a command line argument for the number of attempts.
    error(nargchk(1, 1, nargin))
    
    count = str2num(count);
     
    % Compute pi using monte-carlo simulation
    % pick random points in the unit square. If this falls within a
    % quater circle centred on the origin of radius 1 then count this as
    % a hit. If we do this a number of times then pi/4 = hits/attempts
    hits = 0;
    for i=1:count,
      x = 2*rand-1;  y = 2*rand-1;
      if x^2 + y^2 <= 1, hits = hits + 1; end;
    end;
 
    fprintf('hits = %d\n', hits);
    fprintf('count = %d\n', count);
    fprintf('pi = %1.10f\n', 4*(hits/count));
 
    output = 0;
end
