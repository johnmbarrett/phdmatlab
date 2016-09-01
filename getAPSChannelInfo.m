function [channelInfo,channelIndices] = getAPSChannelInfo
    indices = num2str((1:64)','%02d');
    colIdx = repmat(indices,64,1);
    rowIdx = char(kron(indices,ones(64,1)));
    prefix = repmat('Ch',4096,1);
    infix = repmat('_',4096,1);
    channelInfo = struct('label',cellstr([prefix rowIdx infix colIdx]));
    channelIndices = (1:4096)';
end