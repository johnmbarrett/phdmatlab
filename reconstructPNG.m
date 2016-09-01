function reconstructedData = reconstructPNG(filteredData,filterTypes)
    [nRows,nCols] = size(filteredData);
    
    assert(numel(filterTypes) == nCols);
    
    paddedFilteredData = [zeros(nRows+4,1) [zeros(4,nCols); filteredData]];
    paddedReconstructedData = zeros(size(paddedFilteredData));
        
    for jj = 2:nCols
        filterType = filterTypes(jj);
        
        for ii = 5:nRows
            x = paddedFilteredData(ii,jj);
            a = paddedReconstructedData(ii-4,jj);
            b = paddedReconstructedData(ii,jj-1);
            c = paddedReconstructedData(ii-4,jj-1);
            
            switch(filterType)
                case 0
                    y = x;
                case 1
                    y = mod(x + a,256);
                case 2
                    y = mod(x + b,256);
                case 3
                    y = mod(x + floor(mod((a + b)/2,256)),256); %mod(x + floor((a + b)/2),256);
                case 4
                    y = mod(x + paethPredictor(a,b,c),256);
            end
            
            paddedReconstructedData(ii,jj) = y;
        end
    end
    
    reconstructedData = paddedReconstructedData(5:end,2:end);
end

function pr = paethPredictor(a,b,c)
    p = a+b-c;
    pa = abs(p-a);
    pb = abs(p-b);
    pc = abs(p-c);

    if pa <= pb && pa <= pc
        pr = a;
    elseif pb <= pc
        pr = b;
    else
        pr = c;
    end
end
    