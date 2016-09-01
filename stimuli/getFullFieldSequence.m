function [sequence,nabs,nrel] = getFullFieldSequence(maxLength)
    if nargin < 1
        maxLength = 12;
    end
    
    absIntensities = [0 128 256];
    relIntensities = [-256 -128 128 256];

    nabs = [0 0 0];
    nrel = [0 0 0 0];
    
    function [sequence,nabs,nrel] = getNextIntensity(sequence,nabs,nrel,maxLength)
        if numel(sequence) > maxLength
            error;
        end
        
        if numel(sequence) == maxLength
            if sequence(end) ~= 128
                error;
            end
            
            return;
        end

        ii = sequence(end);
        for jj = setdiff(absIntensities,ii)
            try
                if nabs(absIntensities == jj) + 1 > maxLength/3
                    error;
                end

                if nrel(relIntensities == (jj-ii)) + 1 > maxLength/4
                    error;
                end

                [sequence,nabs,nrel] = getNextIntensity([sequence jj],nabs+(absIntensities == jj),nrel+(relIntensities == (jj-ii)),maxLength);
            catch
                continue;
            end
        end
    end
    
    [sequence,nabs,nrel] = getNextIntensity(128,nabs,nrel,maxLength);
end



