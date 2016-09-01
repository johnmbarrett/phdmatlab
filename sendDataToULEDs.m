function sendDataToULEDs(port,pixels,timings)
    if ~ischar(port) || numel(port) ~= 4 || ~isequal(port(1:3),'COM')
        error 'Invalid serial port specified';
    end
    
    commands = pixelArrayToULEDFrame(pixels);
    
    if ~isnumeric(timings) || ~isvector(timings) || ~all(isfinite(timings)) || any(timings < 1) || numel(timings) ~= size(pixels,3)
        error 'Timings must be a vector of times in milliseconds with length equal to the number of image frames';
    end
    
    com = serial(port,'BaudRate',921600);
    closeCOMPort = onCleanup(@() fclose(com));
    fopen(com);
    
    oldState = pause('on');
    
    for t = 1:numel(timings)
        fprintf(com,'0ffff%s\n',commands(t,:));
        
        pause(timings(t)/1000);
        
        fprintf(com,'9\n');
    end
    
    fclose(com);
    pause(oldState);
end