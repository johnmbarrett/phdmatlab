function mungeAggregateDataForSPSS(outFile,spontData,concData,ns)
    allDatas = { ...
        spontData(:,:,:,1) ...
        spontData(:,:,:,2) ...
        concData(:,:,:,3) ...
        concData(:,:,:,4) ...
        permute(ns(:,:,:,1),[1 3 2]) ...
        };
    
    sheets = {'Firing Rate' 'Oscillation Strength' 'SNR' 'Threshold' 'Responding Cells'};
    
    for ii = 1:5
        mungedData = zeros(0,8);
        
        for jj = 1:3
            data = allDatas{ii}(:,:,jj);
            valid = ~isnan(data(:,1));
            nValid = sum(valid);
            mungedData(end+(1:nValid),2:end) = [repmat(jj,nValid,1) data(valid,:)];
        end
        
        n = size(mungedData,1);
        mungedData(:,1) = 1:n;
        
        xlswrite(outFile,{'Retina' 'Drug' 'Control' '10uM' '20uM' '40uM' '80uM' 'Washout'},sheets{ii},'A1:H1');
        xlswrite(outFile,mungedData,sheets{ii},sprintf('A2:H%d',n+1));
    end
end