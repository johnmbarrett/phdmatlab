function responsiveCells = concatenateResponsiveCells(recordings,varargin)
    allResponsiveCells = zeros(0,2);
    
    options = getopt('prefix='''' concatfun=''union''',varargin{:});
    
    if strcmpi(options.concatfun,'intersect')
        concatfun = @intersect;
    else
        concatfun = @union;
    end
    
    method = 'new'; %#ok<NASGU>
    
    function concatenate(dataFileFormat,convertResponsiveCells)
        for ii = 1:numel(recordings)
            dataFile = sprintf(dataFileFormat,recordings{ii},recordings{ii});

            if ~exist(dataFile,'file')
                continue;
            end

            responsiveCells = load(dataFile,'responsiveCells');
            responsiveCells = responsiveCells.responsiveCells;

            if isempty(responsiveCells)
                continue;
            end

            responsiveCells = convertResponsiveCells(responsiveCells);

            if isempty(allResponsiveCells)
                allResponsiveCells = responsiveCells;
            else
                allResponsiveCells = concatfun(allResponsiveCells,responsiveCells,'rows');
            end
        end
    end
    
    concatenate('%s\\%s_uled_square_responses_newmethod.mat',@(x) x);
        
    if isempty(allResponsiveCells)
        method = 'old'; %#ok<NASGU>
        
        concatenate('%s\\%s_uled_square_responses.mat',@(x) [str2num(x(:,1:2)) str2num(x(:,4))]); %#ok<ST2NM>
    end
    
    responsiveCells = allResponsiveCells;
    
    [~,currentDir] = fileparts(pwd);    
    
    save(sprintf('%s_%s_responsive_cells',currentDir,options.prefix),'responsiveCells','method');
end