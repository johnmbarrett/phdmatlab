function s = arrayToFilenameString(x)
    if isscalar(x)
        s = sprintf('%d',x);
        return;
    end
    
    x = x(:);
    
    if isequal(x,(x(1):x(end))')
        s = sprintf('%d-%d',x(1),x(end));
        return;
    end
    
    s = '';
    
    for ii = 1:numel(x)
        s = sprintf('%s_%d',s,x(ii));
    end
    
    s = s(2:end);
end