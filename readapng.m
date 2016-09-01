function X = readapng(filename)
    fin = fopen(filename);
    
    signature = fread(fin,8);
    
    if ~isequal(signature, [137; 80; 78; 71; 13; 10; 26; 10])
        error 'Not a PNG';
    end
    
    X = [];
    width = NaN;
    height = NaN;
    bpp = NaN;
    colourType = NaN;
    nFrames = NaN;
    colourTypeToColourDimension = [1 NaN 3 NaN 2 NaN 4];
    currentFrame = 0;
%     nextWidth = NaN;
%     nextHeight = NaN;
%     nextXOffset = NaN;
%     nextYOffset = NaN;
%     sequenceNumber = 0;
    
    while true
        lenAndType = fread(fin,8);
        
        if feof(fin)
            break;
        end
        
        len = bytes2num(lenAndType(1:4));
        type = char(lenAndType(5:8));
        
%         disp(type);
        
        if isequal(type','IHDR')
            if len ~= 13
                error 'Not a valid png';
            end
            
            headerData = fread(fin,len);
            
            width = bytes2num(headerData(1:4));
            height = bytes2num(headerData(5:8));
%             nextWidth = width;
%             nextHeight = height;
%             nextXOffset = 0;
%             nextYOffset = 0;
            bpp = headerData(9);
            colourType = headerData(10);
            
            if colourType == 3
                error 'Palette indexed PNGs are not supported';
            end
            
            compressionMethod = headerData(11);
            
            if compressionMethod ~= 0
                error 'Compression methods other than 0 not supported';
            end
            
            filterMethod = headerData(12);
            interlaceMethod = headerData(13);
        elseif isequal(type','acTL')
            if len ~= 8
                error 'Not a valid png';
            end
            
            actlData = fread(fin,len);
            
            nFrames = bytes2num(actlData(1:4));
            
            if isnan(width+height+bpp+colourType+nFrames)
                error 'Not a valid png';
            end
            
            X = zeros(height,width,colourTypeToColourDimension(colourType+1),nFrames);
        % elseif isequal(type','fcTL') % TODO : support fcTL chunks
        elseif isequal(type','IDAT') || isequal(type','fdAT')
            if isequal(type','fDAT')
                sequenceNumber = bytes2num(fread(fin,4));
                len = len - 4;
            else
                sequenceNumber = 1;
            end
            
            frameData = fread(fin,len);
            
            imageData = dunzip(uint8(frameData));
        else
            fseek(fin,len,0);
        end
        
        [~] = fread(fin,4);
    end
    
    X = lens;
end

function num = bytes2num(bytes)
    num = sum(256.^(length(bytes)-1:-1:0)'.*double(bytes));
end
    