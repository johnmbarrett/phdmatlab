function psrhToStimulusResponsePairs(recordings,varargin)
    opts = getopt( ...
            ['latencytype=''firstspike'' ' ...
             'significance=''all'' ' ...
             'rasterfilesuffix=''psrh_abs'' ' ...
             'srpfilesuffix=''srp'' ' ...
             'infovar=2 ' ...
             'confvar=1 ' ...
             'psrhfile=NaN ' ...
             'spiketimesfile=NaN ' ...
             'responsefile=NaN'
            ],varargin{:});
        
    permutation = [1 opts.confvar+1 opts.infovar+1 4];
    
    if ~iscell(recordings)
        recordings = {recordings};
    end
    
    nRecordings = numel(recordings);
    for ee = 1:nRecordings
        recording = recordings{ee};
        %% INTITIALIZATION and defining significant responses
        if ischar(opts.psrhfile)
            psrhFile = opts.psrhfile;
            
            [fileDir,filename] = fileparts(psrhFile);
            
            if isempty(fileDir)
                fileDir = '.';
            end
        else
            [fileDir,filename] = getAnalysisOutputDir(recording);
            
            psrhFile = [fileDir '\' filename '_' opts.rasterfilesuffix '.mat'];
        end
        
        load(psrhFile,'cells','levels','valuess','rasters','stimulusTimings','histograms','edgess');
        cells = [str2num(char(cells(:,1:2))) cells(:,3)]; %#ok<ST2NM>
        
        rasters = permute(rasters,permutation);
        histograms = permute(rasters,permutation); %#ok<NASGU>
        edgess = permute(rasters,permutation); %#ok<NASGU>
        N = numel(rasters{1});

        stimulusTimings = permute(stimulusTimings,[opts.confvar opts.infovar 3]);

        if strcmp(opts.significance,'visual')
            responsiveCells = importdata([fileDir '\lightresponses.txt'],'\t');
        elseif strcmp(opts.significance,'file') && ischar(opts.responsefile) && exist(opts.responsefile,'file')
            responsiveCells = load(opts.responsefile,'responsiveCells');
            responsiveCells = responsiveCells.responsiveCells;
        else
            responsiveCells = cells;
        end
        
        responseIndices = [];
            
        for ii = 1:size(responsiveCells,1)
            index = find(ismember(cells,responsiveCells(ii,:),'rows'));

            if ~isempty(index)
                responseIndices(end+1) = index; %#ok<AGROW>
            end
        end
        
        nResponses = numel(responseIndices);

        assert(nResponses <= size(rasters,1));

        nWidths = levels(opts.confvar);
        nPhases = levels(opts.infovar);
    
        stimuli = nan(nPhases*N,nWidths);
        responses = nan(nPhases*N,2,2,nWidths,nResponses);
        
        %% DATA MUNGING
        widths = valuess{opts.confvar}; %#ok<USENS> %*25/2; %#ok<USENS>
        shifts = valuess{opts.infovar}; %*8; TODO : specify this as a function in varargin
        
        [widths,widthIndices] = sort(widths); %#ok<ASGLU>
        [~,shiftIndices] = sort(shifts);
        
        rasters = rasters(:,widthIndices,shiftIndices,:);
        stimulusTimings = stimulusTimings(widthIndices,shiftIndices,:);

        for ii = 1:nResponses
            cellIndex = responseIndices(ii);
            
            for jj = 1:nWidths
                for kk = 1:nPhases
                    raster = rasters{cellIndex,jj,kk,1};
                    stimulusTiming = stimulusTimings{jj,kk,1};
                    nTrials = numel(raster);
                    
                    stimuli((kk-1)*N+(1:nTrials),jj) = (kk-1)/nPhases; % stimulus is uniformly distributed on the interval [0,1] because it makes my formulae work better
                    
                    for ll = 1:nTrials
                        stimulusStart = stimulusTiming(1,1,ll);
                        row = raster{ll};
                        
                        responseIndex = (kk-1)*N+ll;

                        for mm = 1:size(stimulusTiming,2)
                            validSpikes = row > stimulusTiming(1,mm,ll)-stimulusStart & row <= stimulusTiming(2,mm,ll)-stimulusStart;
                            response = row(validSpikes) - (stimulusTiming(1,mm,ll)-stimulusStart);
                            
                            count = numel(response);
                            
                            if isempty(response)
                                latency = Inf;
                            else
                                latency = response(1);
                            end
                            
                            responses(responseIndex,1,mm,jj,ii) = count;
                            responses(responseIndex,2,mm,jj,ii) = latency;
                        end
                    end
                end
            end
        end
        
        % TODO : new style cells
        save(sprintf('%s\\%s_srp_data',fileDir,filename),'stimuli','responses','cells','responsiveCells','responseIndices','widths','-v7.3');
    end
end