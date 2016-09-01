function [file,filename] = openMCSRawFile(filename)
    fin = fopen(filename);
    
    file = struct;
    file.handle = fin;
%     maybe this wasn't such a good idea..
    file.destructor = onCleanup(@() fclose(file.handle));

    s = fgetl(fin);

    if ~strcmp(s,'MC_DataTool binary conversion')
        warning('Missing header for recording %s, assuming 25kHz sampling rate, signed, all electrodes (row-major layout).  No unit conversion is possible\n', filename);

        file.adcZero = 2^15;
        file.headerSize = 0;
        file.mcdFile = 'Unknown';
        file.nChannels = 60;
        file.sampleRate = 25000;
        file.scaling = 1;
        file.version = 'Unknown';

        electrodeIndices = 1:60;
        electrodeChannels = [repmat('El_',60,1) num2str([21:10:71 repmat(12:10:82,1,6)+kron(0:5,ones(1,8)) 28:10:78]')];
    else
        file.version = fscanf(fin,'Version %s\n');

        mcdFile = fgetl(fin);

        [~,~,~,~,tokens] = regexp(mcdFile,'^MC_REC file = "(.*)"$');

        file.mcdFile = tokens{1}{1};

        file.sampleRate = fscanf(fin,'Sample rate = %d\n');
        file.adcZero = fscanf(fin,'ADC zero = %d\n');
        file.scaling = fscanf(fin,'El = %fµV/AD\n');

        streams = fscanf(fin,'Streams = %s\n');
        streams = regexp(streams,';','split');
        
        file.nStreams = numel(streams);

        % TODO : other types of streams
        electrodeIndices = find(strncmp('El',streams,2));
        electrodeChannels = streams(electrodeIndices);
        electrodeChannels = vertcat(electrodeChannels{:});
        electrodeChannels = electrodeChannels(:,end-1:end);

        s = fgetl(fin);

        assert(strcmp(s,'EOH'),'Invalid header format for recording %s',filename);

        file.headerSize = ftell(fin);
    end

    file.electrodes = struct( ...
        'index',    num2cell(electrodeIndices(:)), ...
        'label',    cellstr(electrodeChannels) ...
        );
    
    fseek(fin,0,1);
    
    file.dataSize = ftell(fin)-file.headerSize;
    file.nSamples = ceil(file.dataSize/(2*file.nStreams)); % /2 because each sample is two bytes
    
    fseek(fin,file.headerSize,-1);
end