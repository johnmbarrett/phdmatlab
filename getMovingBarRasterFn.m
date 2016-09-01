function f = getMovingBarRasterFn(rfFileIndex)
    f = @(ax,yy,recording,channel,cluster,subscript,valuess) movingBarRasterFn(ax,yy,recording,channel,cluster,subscript,valuess,rfFileIndex);
end

function movingBarRasterFn(~,yy,~,channel,cluster,subscript,valuess,rfFileIndex)
    recordings = initRecordings;
    rfRecording = recordings(rfFileIndex);
    rfDir = getAnalysisOutputDir(rfRecording);
    rfFile = [rfDir '\rf_' rfRecording.dataFile '_channel_' channel '_cluster_' cluster '_clustered_phototrigger_all_spikes.mat'];
    
    try
        load(rfFile,'centreX','centreY','sd');
    catch err
        logMatlabError(err);
        return;
    end
    
    starts = valuess{1};
    centre = mean(starts);
    
    start = starts(subscript{2},:);
    finish = start+2*(centre-start);
    speed = valuess{2}(subscript{3});
    width = valuess{4}(subscript{5});
    
    [Y,X] = ndgrid(1:500,1:500);
    
    theta = atan2(finish(2)-start(2),finish(1)-start(1));
    
    if ismember(theta,[-pi 0 pi]);
        slope = Inf;
    else
        slope = tan(theta + pi/2);
    end
        
    trajectoryLength = ceil(sqrt(sum((finish-start).^2)));
    
    tEnter = NaN;
    tExit = NaN;
    for t = 0:trajectoryLength
        x0 = start(1)+(t)*cos(theta);
        y0 = start(2)+(t)*sin(theta);
        
        if ~isfinite(slope)
            edge = X == round(x0);
        else
            edge = abs(slope*(X-x0)+y0-Y) < 2;
        end
        
        inRF = ((X-centreX).^2 + (Y-centreY).^2) <= (2*sd)^2;
        
        if isnan(tEnter) && any(any(edge & inRF))
            tEnter = t-ceil(width/2);
        end
        
        if ~isnan(tEnter) && isnan(tExit) && ~any(any(edge & inRF))
            tExit = t+ceil(width/2);
            break;
        end
    end
    
    if isnan(tEnter)
        return;
    end
    
    line(repmat(tEnter/(60*speed),1,2),yy,'LineStyle','--','Color',[0.5 0.5 0.5]);
    
    if ~isnan(tExit)
        line(repmat(tExit/(60*speed),1,2),yy,'LineStyle','--','Color',[0.5 0.5 0.5]);
    end
end