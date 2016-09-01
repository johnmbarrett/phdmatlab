function logMatlabError(err)
    fprintf('Error %s: %s\n',err.identifier,err.message);
        
    for ii = 1:numel(err.stack)
        fprintf('On line %d of file %s\n',err.stack(ii).line,err.stack(ii).file);
    end
end