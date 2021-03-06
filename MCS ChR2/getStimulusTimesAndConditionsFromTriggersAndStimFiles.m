function [onsetss,offsetss,sequences] = getStimulusTimesAndConditionsFromTriggersAndStimFiles(triggerFile,stimSpecs,specIndices,nTriggerss,types,outputFilenames,stimTimesFilePrefix)
    nSpecs = numel(stimSpecs);
    
    if nargin < 7
        stimTimesFilePrefix = '';
    else
        stimTimesFilePrefix = [stimTimesFilePrefix '_'];
    end
    
    if nargin < 6
        outputFilenames = arrayfun(@(n) sprintf('sequence %d.mat',n),1:nSpecs,'UniformOutput',false);
    end
    
    if nargin < 5
        types = cell(nSpecs,1);
        
        for ii = 1:nSpecs
            stimulusTypes = unique({stimSpecs{ii}.type});
            stimulusTypes = setdiff(stimulusTypes,{'b'});
            
            if numel(stimulusTypes) > 1
                error('Multiple stimulus types used in sequence %d.  Please specify which stimulus type you wish to include in the sequence.\n',ii);
            end
            
            if numel(stimulusTypes) < 1
                warning('Sequence %d only contains blank stimuli, why are you even analysing it?\n',ii)
                types{ii} = 'b';
                continue;
            end
            
            types(ii) = stimulusTypes;
        end
    end
        
    if isnumeric(triggerFile) && isvector(triggerFile)
        triggerTimes = triggerFile;
    else
        load(triggerFile,'vsyncTimes','recordingStartTime');
        triggerTimes = vsyncTimes-recordingStartTime;
    end
    
    assert(numel(triggerTimes) == sum(nTriggerss),[ ...
        'Mismatch between actual and parsed number of triggers.  Possible explanations:\n' ...
        ' * You chose the wrong trigger file\n' ...
        ' * You chose the wrong set of stim files\n' ...
        ' * Your trigger signal was corrupt\n' ...
        ' * You started the recording late\n' ...
        ' * You stopped the recording early\n']);
    
    assert(numel(specIndices) == numel(nTriggerss),'Each stim file should have both a stim spec and a number of triggers');
    
    cTriggers = [0; cumsum(nTriggerss)];
    
    nRecordings = numel(nTriggerss);
    
    onsetss = cell(nRecordings,1);
    offsetss = cell(nRecordings,1);
    
    for ii = 1:nRecordings
        specIndex = specIndices(ii);
        stimSpec = stimSpecs{specIndex};
        type = types{specIndex};
        
        relevantTriggers = triggerTimes(cTriggers(ii)+1:cTriggers(ii+1));
        
        onsetIndices = vertcat(stimSpec(strcmp(type,{stimSpec.type})).index);
        onsetss{ii} = relevantTriggers(onsetIndices);
        
        offsetIndices = vertcat(stimSpec(strcmp('b',{stimSpec.type})).index);
        
        if isempty(offsetIndices)
            offsetss{ii} = [onsetss{ii}(2:end); onsetss{ii}(end)+stimSpec(end).duration/1000];
        else
            offsetss{ii} = relevantTriggers(offsetIndices);
        end
    end
    
    save([stimTimesFilePrefix 'stimulusTimes.mat'],'onsetss','offsetss');
    
    sequences = cell(nSpecs,1);
    
    for ii = 1:nSpecs
        type = types{ii};
        
        stimSpec = stimSpecs{ii};
        stimSpec = stimSpec(strcmp(type,{stimSpec.type}));
        
        nStimuli = numel(stimSpec);
        
        assert(all(cellfun(@numel,onsetss(specIndices == ii)) == nStimuli),'Mismatch between number of stimuli specific and number of stimulus onsets found');
        assert(all(cellfun(@numel,offsetss(specIndices == ii)) == nStimuli),'Mismatch between number of stimuli specific and number of stimulus offsets found');
        
        switch type
            case 'b'
                error('This function has not been implemented yet.  Yell at John to implement it, or try writing your own based on the examples below!');
            case 'c'
                error('There are no contracting rings.  Only rings that expand inwards.');
            case 'e'
                error('This function has not been implemented yet.  Yell at John to implement it, or try writing your own based on the examples below!');
            case 'g'
                error('This function has not been implemented yet.  Yell at John to implement it, or try writing your own based on the examples below!');
            case 'i'
                error('This function has not been implemented yet.  Yell at John to implement it, or try writing your own based on the examples below!');
            case 'm'
                varNames = {'direction' 'period' 'width'};
            case 'r'
                varNames = {'X' 'Y' 'width' 'height' 'duration'};
            case 'w'
                error('This function has not been implemented yet.  Yell at John to implement it, or try writing your own based on the examples below!');
        end
        
        nVars = numel(varNames);
        
        sequence = struct('conditions',[],'conditionOrder',[],'varNames',[]);
        conditionOrder = zeros(nStimuli,nVars);
        
        for jj = 1:nVars
            varName = varNames{jj};
            
            if any(strcmp(varName,{'direction' 'path'}))
                values = {stimSpec.(varName)};
                values = values(:);
            else
                values = vertcat(stimSpec.(varName));
            end
            
            [sequence.(varNames{jj}),~,conditionOrder(:,jj)] = unique(values);
        end
        
        degenerate = all(conditionOrder == 1);
        
        goodVars = find(~degenerate);
        nGoodVars = numel(goodVars);
        
        for jj = nGoodVars:-1:2
            for kk = jj-1:-1:1
                if isequal(conditionOrder(:,goodVars(jj)),conditionOrder(:,goodVars(kk)))
                    degenerate(jj) = true;
                    break;
                end
            end
        end
        
        sequence.varNames = varNames(~degenerate);
        [sequence.conditions,~,sequence.conditionOrder] = unique(conditionOrder(:,~degenerate),'rows');
        
        assert(isequal(conditionOrder(:,~degenerate),sequence.conditions(sequence.conditionOrder,:)),'Error converting stim spec into condition list');
        
        save(outputFilenames{ii},'-struct','sequence',varNames{:},'conditions','conditionOrder','varNames');
        
        sequences{ii} = sequence;
    end
end