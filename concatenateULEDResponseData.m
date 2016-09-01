function concatenateULEDResponseData(recordings,concentrations,voltages,varargin)
    if numel(recordings) ~= numel(concentrations) || numel(concentrations) ~= numel(voltages)
        error('Number of recordings, drug concentrations and uLED voltages must match');
    end
    
    recordings = recordings(:);
    
    options = getopt('format=''graphpad''',varargin{:});
    
    allResponsiveCells = zeros(0,4);
    allThresholds = zeros(0,5);
    allNSpikes = zeros(0,6);
    allAmplitudes = zeros(0,6);
    
    allCellIndices = zeros(87,10,2);
    
    nRecordings = numel(recordings);
    for ii = 1:nRecordings
        recording = recordings{ii};
        [filedir,filename] = getAnalysisOutputDir(recording);
        
        responseFile = sprintf('%s\\%s_uled_square_responses',filedir,filename);
        load(responseFile);
        
        c = concentrations(ii);
        v = voltages(ii);
        
        nCells = size(responsiveCells,1);
        allResponsiveCells = [allResponsiveCells; nCells c v ii]; %#ok<*AGROW>
        
        if nCells == 0
            continue;
        end
        
        cellIndices = zeros(nCells,1);
        
        for jj = 1:nCells
            cellID = sscanf(responsiveCells(jj,:),'%d %d',[1 2]);
            channel = cellID(1);
            cluster = cellID(2);
            voltage = (v == 6) + 1;
            
            if allCellIndices(channel,cluster,voltage) > 0
                cellIndices(jj) = allCellIndices(channel,cluster,voltage);
            else
                index = max(max(allCellIndices(:,:,voltage)))+1;
                cellIndices(jj) = index;
                allCellIndices(channel,cluster,voltage) = index;
            end
        end
        
        n = numel(thresholdss); %#ok<*NODEF>
        allThresholds = [allThresholds; thresholdss(:) cellIndices c*ones(n,1) v*ones(n,1) ii*ones(n,1)];
        
        if ~exist('pws','var');
            pws = [5 10 25 50 75 100];
        end
        
        nPWs = numel(pws);
        nTrials = size(nSpikess,1)/nCells;
        
        nSpikess = squeeze(median(reshape(nSpikess,nTrials,nCells,nPWs)));
        
        n = nCells*nPWs;
        nSpikess = reshape(nSpikess,n,1);
        
        allCells = repmat(cellIndices,nPWs,1);
        allPWs = kron(pws(:),ones(nCells,1));
        allNSpikes = [allNSpikes; nSpikess allCells allPWs c*ones(n,1) v*ones(n,1) ii*ones(n,1)];

        amplitudess = squeeze(median(reshape(amplitudess,nTrials,nCells,nPWs)));
        amplitudess = reshape(amplitudess,n,1);
        
        allAmplitudes = [allAmplitudes; amplitudess allCells allPWs c*ones(n,1) v*ones(n,1) ii*ones(n,1)];
    end
    
    [~,currentDir] = fileparts(pwd);
    
    dataFile = sprintf('%s_uled_square_responses_%s.xlsx',currentDir,options.format);
    
    if exist(dataFile,'file');
        delete(dataFile);
    end
    
    if strcmpi(options.format,'column')
        xlswrite(dataFile,{'# Cells' 'Concentration' 'Voltage' 'Recording'},'# Responsive Cells','A1:D1');
        xlswrite(dataFile,allResponsiveCells,'# Responsive Cells',sprintf('A2:D%d',size(allResponsiveCells,1)+1));

        xlswrite(dataFile,{'Threshold' 'Cell' 'Concentration' 'Voltage' 'Recording'},'Thresholds','A1:D1');
        xlswrite(dataFile,allThresholds,'Thresholds',sprintf('A2:E%d',size(allThresholds,1)+1));

        xlswrite(dataFile,{'# Spikes' 'Cell' 'Pulse Width' 'Concentration' 'Voltage' 'Recording'},'Spikes Per Pulse','A1:E1');
        xlswrite(dataFile,allNSpikes,'Spikes Per Pulse',sprintf('A2:F%d',size(allNSpikes,1)+1));

        xlswrite(dataFile,{'Amplitude' 'Cell' 'Pulse Width' 'Concentration' 'Voltage' 'Recording'},'Spike Amplitudes','A1:E1');
        xlswrite(dataFile,allAmplitudes,'Spike Amplitudes',sprintf('A2:F%d',size(allAmplitudes,1)+1));
    elseif strcmpi(options.format,'graphpad');
        xlswrite(dataFile,{'Concentration' '5V' '6V'},'# Responsive Cells','A1:C1');
        
        uniqueConcs = unique(concentrations)';
        n = numel(uniqueConcs);
        
        try
            xlswrite(dataFile,[... 
                uniqueConcs ...
                allResponsiveCells(allResponsiveCells(:,3) == 5 & allResponsiveCells(:,4) < nRecordings-1,1) ...
                allResponsiveCells(allResponsiveCells(:,3) == 6 & allResponsiveCells(:,4) < nRecordings-1,1) ...
                ],'# Responsive Cells',sprintf('A2:C%d',n+1));
        catch err %#ok<NASGU>
            warning('Not all voltages were recorded at all concentrations');
        end
        
        n6VRecordings = sum(voltages == 6);
        lastColumn = char('A'+n6VRecordings-1);
        
        i6VRecordings = find(voltages == 6);
        xlswrite(dataFile,recordings(i6VRecordings(:))','Thresholds',sprintf('A1:%s1',lastColumn));
        
        maxNCells = max(max(allCellIndices(:,:,2)));
        graphPadThresholds = nan(maxNCells,n6VRecordings);
        
        for ii = 1:n6VRecordings
            rows = allThresholds(:,5) == i6VRecordings(ii);
            graphPadThresholds(allThresholds(rows,2),ii) = allThresholds(rows,1);
        end
        
        xlswrite(dataFile,graphPadThresholds,'Thresholds',sprintf('A2:%s%d',lastColumn,maxNCells+1));
        
        header = [{'Pulse Width'} recordings(i6VRecordings(:))'];
        
        pws = unique(allPWs);
        nPWs = numel(pws);
        
        allSheets = {'Spikes Per Pulse (all)', 'Spike Amplitude (all)'};
        maxPWSheets = {'Spikes Per Pulse (max PW)', 'Spike Amplitude (maxPW)'};
        
        datas = zeros([size(allNSpikes) 2]);
        datas(:,:,1) = allNSpikes;
        datas(:,:,2) = allAmplitudes;
        
        for hh = 1:2
            xlswrite(dataFile,header(1),allSheets{hh},'A1:A1');
            xlswrite(dataFile,header(2:end),maxPWSheets{hh},sprintf('A1:%s1',lastColumn));
            xlswrite(dataFile,pws,allSheets{hh},sprintf('A2:A%d',1+nPWs));

            graphPadDataAll = nan(nPWs,maxNCells*n6VRecordings);
            graphPadDataMaxPW = nan(maxNCells,n6VRecordings);

            for ii = 1:n6VRecordings
                rows = datas(:,6,hh) == i6VRecordings(ii);
                data = datas(rows,1,hh);
                nCells = size(data,1)/nPWs;
                data = reshape(data,nCells,nPWs)';
                cellIndices = datas(rows,2,hh);
                assert(isequal(cellIndices,repmat(cellIndices(1:nCells),nPWs,1)));
                cellIndices = cellIndices(1:nCells);
                graphPadDataAll(:,maxNCells*(ii-1)+(cellIndices)) = data;
                graphPadDataMaxPW(cellIndices,ii) = data(end,:)';
                column = getExcelColumn((ii-1)*maxNCells+1);
                xlswrite(dataFile,header(ii+1),allSheets{hh},sprintf('%s1:%s1',column,column));
            end

            xlswrite(dataFile,graphPadDataAll,allSheets{hh},sprintf('B2:%s%d',getExcelColumn(n6VRecordings*maxNCells),1+nPWs));
            xlswrite(dataFile,graphPadDataMaxPW,maxPWSheets{hh},sprintf('A2:%s%d',lastColumn,1+maxNCells));
        end
    end
end  