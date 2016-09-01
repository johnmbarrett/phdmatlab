function [figs,lines,edges,hists,trials,subs,nSpikess] = rasterPlot(figs,varargin) %spikeTimes,stimulusTimes,conditions,tmin,tmax,isPSTH,bw,groups,varNames,isHz)
    isHandle = all(ishandle(figs));
    
    noPlot = false;
    if isHandle
        args = varargin;
        handleTypes = get(figs,'type');
        
        if ~all(strcmp('figure',handleTypes))
            error('Unknown handle type');
        end
        
        nFigs = numel(figs);
        
        for ii = 1:numel(figs)
            set(figs(ii),'Visible','off');
            clf(figs(ii));
        end
    elseif (islogical(figs) && ~any(figs)) || any(isnan(figs))
        noPlot = true;
        args = varargin;
        nFigs = 0;
    else
        args = {figs varargin{:}}; %#ok<CCAT>
        figs = figure('Visible','off');
        nFigs = 1;
    end
    
    nargin = numel(args);
    
    if nargin < 2
        error('rasterPlot needs, at bare minimum, vectors of spike times & stimulus times');
    end
    
    if nargin < 10
        isHz = true;
    else
        isHz = all(logical(args{10}));
    end
    
    spikeTimes = args{1};
    stimulusTimes = args{2};
    
    if nargin < 9
        varNames = NaN;
    else
        varNames = args{9};
    end
    
    if nargin < 8 || isempty(args{8})
        groups = ones(size(stimulusTimes));
    else
        groups = args{8};
    end
    
    nVars = size(groups,2)+1-isequal(groups,ones(size(groups)));
    if ~iscell(varNames) || numel(varNames) ~= nVars
        varNames = cellstr([repmat('Var ',nVars) num2str((1:nVars)')]);
    end
    
    uniqueGroups = unique(groups,'rows');
    nGroups = size(uniqueGroups,1);
    
    if noPlot
        figs = NaN;
    else
        figs((nFigs+1):nGroups) = zeros(nGroups-nFigs,1);

        for ii = (nFigs+1):nGroups
            figs(ii) = figure('Visible','off');
        end
    end
        
    if nargin < 7
        bw = 0.01;
    else
        bw = args{7};
    end

    if nargin < 6
        isPSTH = false;
    else
        isPSTH = args{6};
    end

    if nargin < 5
        if nargin < 4
            ts = [0; 1];
        else
            ts = args{4};
        end
    else
        ts = [args{4}; args{5}];
        ts = ts(isfinite(ts));
    end
    
    if ~ismember(0,ts)
        ts(end+1) = 0;
    end
    
    ts = sort(ts);
    ts = ts(:);
    
    if nargin < 3
        conditions = ones(size(stimulusTimes));
    else
        conditions = args{3};
    
        if size(conditions,1) ~= size(stimulusTimes,1)
            error('Number of condition specifiers must match the number of stimuli');
        end
        
        if iscell(varNames{1}) && size(conditions,2) ~= numel(varNames{1})
            error('Number of sub-conditions must match number of variable names');
        end 
    end

    % suppress annoying warning
    if iscell(conditions)
        uniqueConditions = unique(conditions);
    else
        uniqueConditions = unique(conditions,'rows');
    end
    
    nConditions = size(uniqueConditions,1);
    
    trials = zeros(nConditions,nGroups);
    
    [rows,cols] = subplots(nConditions);
    
    lines = cell(nConditions,nGroups);
    nSpikess = cell(nConditions,nGroups);
    
    if isPSTH
        edges = ((floor(ts(1)/bw)*bw:bw:bw*ceil(ts(end)/bw)))';
        nEdges = numel(edges);
        hists = zeros(nEdges-1,nConditions,nGroups);
    else
        edges = NaN;
        hists = NaN;
    end
    
    for ii = 1:numel(stimulusTimes);
        if iscell(varNames{1})
            condition = find(ismember(uniqueConditions,conditions(ii,:),'rows'));
        elseif iscell(uniqueConditions)
            condition = find(strcmp(conditions(ii),uniqueConditions));
        else
            condition = find(uniqueConditions == conditions(ii));
        end
        
        group = ismember(uniqueGroups,groups(ii,:),'rows');
        
%         figure(figs(group));
%         
%         subplot(rows,cols,condition);
        
        trials(condition,group) = trials(condition,group)+1;
        
        stimulusTime = stimulusTimes(ii);
        spikesInTrial = spikeTimes(spikeTimes >= stimulusTime+ts(1) & spikeTimes < stimulusTime+ts(end))-stimulusTime;
        
        lines{condition,group}{end+1} = spikesInTrial;
        
        if isPSTH
            if isempty(spikesInTrial)
                nSpikes = zeros(1,nEdges-1);
            else
                nSpikes = reshape(histc(spikesInTrial,edges),1,nEdges);
                nSpikes = nSpikes(1:nEdges-1); % last bin is always empty
            end
            
            hists(:,condition,group) = hists(:,condition,group) + nSpikes';
        else
            nSpikes = numel(spikesInTrial);
        end
        
        if isempty(nSpikess{condition,group})
            nSpikess{condition,group} = nSpikes;
        else
            nSpikess{condition,group}(end+1,:) = nSpikes;
        end
    end
    
    if isPSTH
        if isHz
            hists = hists./(repmat(reshape(trials,[1 nConditions nGroups]),nEdges-1,1)*bw);
        end
        
        maxH = max(max(max(hists)));
    else
        maxH = 0;
    end
    
    if noPlot
        subs = NaN;
        return;
    end
    
    subs = zeros(nConditions,nGroups);
    
    for hh = 1:nGroups
        set(0,'CurrentFigure',figs(hh));
        
        for ii = 1:nConditions
            if trials(ii,hh) == 0
                continue;
            end
            
            subs(ii,hh) = subplot(rows,cols,ii);
            xlim(ts([1 end]));

            trialOffset = maxH;

            if isPSTH && maxH > 0
                trialHeight = maxH/trials(ii,hh);
                yy = [0 2*maxH];
            else
                trialHeight = 1;
                yy = [0 trials(ii,hh)];
            end
            
            ylim(yy);

            hold on;

            if numel(ts > 2)
                line(repmat(ts(2:end-1)',2,1),repmat(yy',1,numel(ts)-2),'Color','k','LineStyle','--');
            end

            for jj = 1:trials(ii,hh)
                spikesInTrial = lines{ii,hh}{jj};

                if isempty(spikesInTrial)
                    continue;
                end

                line( ...
                    repmat(spikesInTrial(:)',2,1), ...
                    repmat((jj-[1; 0])*trialHeight+trialOffset,1,numel(spikesInTrial)), ...
                    'Color',    'b' ...
                    );
            end

            if isPSTH
                bar(edges(1:end-1),hists(:,ii,hh),'histc');
            end
            
            if ~iscell(varNames{1})
                conditionVars = varNames(1);
            else
                conditionVars = varNames{1};
            end
            
            tit = '';
            for jj = 1:numel(conditionVars)
                if iscell(uniqueConditions)
                    tit = sprintf('%s, %s = %s',tit,conditionVars{jj},uniqueConditions{ii,jj});
                elseif mod(uniqueConditions(ii,jj),1) == 0
                    tit = sprintf('%s, %s = %d',tit,conditionVars{jj},uniqueConditions(ii,jj));
                else
                    tit = sprintf('%s, %s = %f',tit,conditionVars{jj},uniqueConditions(ii,jj));
                end
            end
            title(tit(3:end));
            
            if ii == (rows-1)*cols+1
                xlabel('Time/s');
                
                if isPSTH
                    if isHz
                        ylabel('Firing Rate/Hz');
                    else
                        ylabel('Spikes/Bin');
                    end
                else
                    ylabel('Trials');
                end
            end
        end
        
        supertitle = '';
        
        for ii = 2:nVars
            if mod(uniqueGroups(hh,ii-1),1) == 0
                dataType = 'd';
            else
                dataType = 'f';
            end
            
            supertitle = sprintf(['%s' repmat(', ',1,ii > 2) '%s = %' dataType],supertitle,varNames{ii},uniqueGroups(hh,ii-1));
        end
        
        suptitle(supertitle);
    end
end