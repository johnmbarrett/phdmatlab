function X = makeSTAToolkitStruct(s,R,trialLength,label,samplingFrequency,startTime)
    if nargin < 6
        startTime = 0;
    end
    
    [us,~,iu] = unique(s);
    M = numel(us);
    X = struct('M',int32(M),'N',int32(1));
    X.sites = struct('label',{{label}},'recording_tag',{{'episodic'}},'time_scale',1,'time_resolution',1/samplingFrequency,'si_unit','none','si_prefix',1);
    
    P = int32(accumarray(iu,ones(size(s))));
    
    X.categories = struct('label',arrayfun(@(u) {u},us,'UniformOutput',false),'P',num2cell(P));
    
    for kk = 1:M
        r = R(iu == kk);
        Q = cellfun(@(x) int32(numel(x)),r,'UniformOutput',false);
        X.categories(kk).trials = struct('start_time',startTime,'end_time',trialLength,'Q',Q,'list',r);
    end
end