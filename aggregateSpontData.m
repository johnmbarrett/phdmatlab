function data = aggregateSpontData(experiments,drugs)
    nDrugs = numel(drugs);
    
    data = arrayfun(@(x) zeros(0,6,3),1:nDrugs,'UniformOutput',false);
    
    suffixes = {'firing_rates' 'fft_data' 'fft_data'};
    sheets = {'Sheet1' 'Peak Power (local peak)' 'Peak Frequency'};
    
    for ii = 1:numel(experiments)
        expDir = sprintf('JBOG%04d',experiments(ii));
        
        for jj = 1:numel(drugs)
            todo = 3;
            done = 0;
            nn = 1;
            while done < todo;
                for kk = 1:3
                    dataFile = sprintf('%s\\%s_%s_%s.xlsx',expDir,expDir,drugs{jj},suffixes{kk});

                    if ~exist(dataFile,'file')
                        dataFile = dataFile(1:end-1); % sometimes it's an xls not xlsx?

                        if ~exist(dataFile,'file')
                            dataFile = sprintf('%s\\%s_%s_%d_%s.xlsx',expDir,expDir,drugs{jj},nn,suffixes{kk}); % handle the double flupirtine experiment
                            
                            if ~exist(dataFile,'file')
                                dataFile = dataFile(1:end-1); % sometimes it's an xls not xlsx?
                                
                                if ~exist(dataFile,'file')
                                    todo = todo - 1;
                                    continue;
                                end
                            end
                            
                            todo = todo + 1;
                        end
                    end

                    X = xlsread(dataFile,sheets{kk},'A2:F61');

                    if kk == 1
                        validChannels = X(:,1) ~= 0;
                    end

                    X = X(validChannels,:);
                    
                    if kk < 3
                        X = X./repmat(X(:,1),1,6);
                    end

                    data{jj}(end+(kk == 1),:,kk) = median(X);
                    
                    done = done + 1;
                end
                
                nn = nn + 1;
            end
        end
    end
end