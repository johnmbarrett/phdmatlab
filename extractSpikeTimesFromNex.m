function extractSpikeTimesFromNex(infiles,outfiles,isSorted)
    if nargin < 3
        isSorted = true;
    end
    
    nonExistent = @(file) ~exist(file,'file');
    if nargin < 1 || isempty(infiles) || (ischar(infiles) && nonExistent(infiles)) || (iscell(infiles) && any(cellfun(nonExistent,infiles)))
        infiles = uipickfiles('FilterSpec','*.nex');
    elseif ischar(infiles)
        infiles = {infiles};
    elseif ~iscell(infiles)
        error('First argument must be empty, a path to a single input file or a cell array of paths to input files');
    end
    
    nFiles = numel(infiles);
    
    if nargin < 2 || isempty(outfiles)
        if nFiles == 1
            [outfile,outpath] = uiputfile('*.mat');
            outfiles = {[outpath outfile]};
        else
            outfiles = {cellfun(@(s) strrep(s,'.nex',' spike times.mat'))};
        end
    elseif ischar(outfiles)
        outfiles = {outfiles};
    elseif ~iscell(outfiles)
        error('Second argument must be empty, a path to a single output file (if a single input file was specified) or a cell array of paths to output files for each input file');
    end
    
    assert(numel(outfiles) == nFiles,'There must be one output file for every input file');
        
    for hh = 1:nFiles
        infile = infiles{hh};
        
        nex = readNexFile(infile);
        unitNames = cellfun(@(wave) wave.name,nex.waves,'UniformOutput',false);
        spikeTimess = cellfun(@(wave) wave.timestamps,nex.waves,'UniformOutput',false);
        
        if isSorted
            % need this step if exported templates as single-unit waveforms
            notTemplates = ~cellfun('isempty',strfind(unitNames,'wf'));
            unitNames = unitNames(notTemplates); %#ok<NASGU>
            spikeTimess = spikeTimess(notTemplates); %#ok<NASGU>
        end
        
        save(outfiles{hh},'unitNames','spikeTimess');
    end
end