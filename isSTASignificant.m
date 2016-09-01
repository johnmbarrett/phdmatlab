function [significant, channels, clusters, zmax, zmin] = isSTASignificant(recording,varargin)
    fileDir = getAnalysisOutputDir(recording);
    
    options = getopt('srrf=false',varargin{:});
    
    if options.srrf
        srrfs = dir([fileDir '\srrf*.mat']);
        srrfs = {srrfs.name}';
        sig = false(size(srrfs));
        channels = zeros(size(srrfs));
        clusters = zeros(size(srrfs));
        for ii = 1:numel(srrfs)
        srrf = srrfs{ii};
        load([fileDir '\' srrf]);
        sig(ii) = significant; %#ok<NODEF>
        [~,~,~,~,channel] = regexp(srrf,'channel_([0-9]+)','once');
        channels(ii) = str2double(channel{1});
        [~,~,~,~,cluster] = regexp(srrf,'cluster_([0-9]+)','once');
        clusters(ii) = str2double(cluster{1});
        end
        chcl = [channels clusters];
        visFile = [fileDir '\' recording.dataFile '_stavis.xlsx'];
        if ~exist(visFile,'file')
            disp(sum(sig));
            disp([channels(sig) clusters(sig)]);
            return;
        end
        vis = xlsread(visFile,['A2:E' num2str(numel(srrfs)+1)]);
        vis(isnan(vis)) = 0;
        isequal(chcl,vis(:,1:2))
        pos = vis(:,3);
        neu = vis(:,4);
        neg = vis(:,5);
        tp = sum(sig & pos) %#ok<NOPRT>
        fp = sum(sig & neg) %#ok<NOPRT>
        tn = sum(~sig & neg) %#ok<NOPRT>
        fn = sum(~sig & pos) %#ok<NOPRT>
        lib = sum(sig & neu) %#ok<NOPRT>
        con = sum(~sig & neu) %#ok<NOPRT>
        tp + fp + fn + tn + lib + con %#ok<NOPRT>
        [channels(~sig & pos) clusters(~sig & pos)] %#ok<NOPRT>
        [channels(sig & neg) clusters(sig & neg)] %#ok<NOPRT>
        return;
    end
    
    gifs = dir([fileDir '\sta_*.gif']);
    gifs = {gifs.name}';
    
    nGifs = numel(gifs);
    channels = zeros(nGifs,1);
    clusters = zeros(nGifs,1);
    zmax = zeros(nGifs,1);
    zmin = zeros(nGifs,1);
%     alphas = zeros(nGifs,1);  
    alpha = 0.05; %2*normcdf(-4.5,0,1); 
    significant = false(nGifs,1);
    maxSigP = zeros(nGifs,1);
    
    for ii = 1:nGifs
        gif = gifs{ii};
        
        pixels = double(squeeze(imread([fileDir '\' gif],'frames','all')));
        N = numel(pixels);
        Z = reshape(pixels,N,1);
        Z = (Z - mean(Z))/std(Z);
        zmax(ii) = max(Z);
        zmin(ii) = min(Z);
        
        P = 2*normcdf(-abs(Z));
%         significant(ii) = any(P < alpha);
        
        [hypotheses,maxSigP(ii)] = fdrcorrect(P,alpha);
        significant(ii) = any(hypotheses); %sum(hs)/numel(hs) >= alpha;
        
%         alphas(ii) = norminv(0.025/N);
        
        [~,~,~,~,channel] = regexp(gif,'channel_([0-9]+)','once');
        channels(ii) = str2double(channel{1});
        
        [~,~,~,~,cluster] = regexp(gif,'cluster_([0-9]+)','once');
        clusters(ii) = str2double(cluster{1});
    end
    
%     significant = zmax > -alphas | zmin < alphas;
    
    save([fileDir '\' recording.dataFile '_stastats.mat'],'significant','channels','clusters','zmax','zmin'); %,'maxSigP');
end