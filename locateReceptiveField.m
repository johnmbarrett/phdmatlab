function locateReceptiveField(recording,varargin)
    channelInfo = getMetaInfo(recording.index,'spont',true);
    channelIndices = parseChannelOptions(recording.index,varargin{:});
    fileDir = getAnalysisOutputDir(recording);
    
    for ii = 1:numel(channelIndices)
        channelIndex = channelIndices(ii);
        channelLabel = channelInfo(channelIndex).label;
        
        stas = dir([fileDir '\sta_ ' recording.dataFile '_channel_' channelLabel '_cluster_*.gif']);
        stas = {stas.name};
        
        for jj = 1:numel(stas)
            sta = [fileDir '\' stas{jj}];
            
            I = squeeze(double(imread(sta,'frames','all')));
            J = reshape(I,numel(I),1);
            
            imax = find(I == max(J));
            [y,x,t] = ind2sub(size(I),imax);
            centreY = ceil(mean(y));
            centreX = ceil(mean(x));
            centreT = ceil(mean(t));
            
            Z = (I-mean(J))/std(J);
           
            sds = (1:(1/25):max(size(I,1),size(I,2)))';
            mses = zeros(size(sds));
            maxZ = max(max(max(Z)));
            
            [Y,X] = ndgrid(1:96,1:96);
            Y = reshape(Y,numel(Y),1);
            X = reshape(X,numel(X),1);
            
            for kk = 1:numel(sds)
                tic;
                sd = sds(kk);
                rf = maxZ*mvnpdf([Y X],[centreY centreX],[sd sd])/mvnpdf([centreY centreX],[centreY centreX],[sd sd]);
                mses(kk) = mean((reshape(squeeze(Z(:,:,centreT)),size(rf))-rf).^2);
                toc;
            end
            
            sd = sds(mses == min(mses));
            rf = maxZ*mvnpdf([Y X],[centreY centreX],[sd sd])/mvnpdf([centreY centreX],[centreY centreX],[sd sd]);
            rf = reshape(rf,size(I,1),size(I,2)); %#ok<NASGU>
            centreX = centreX+recording.textureRect(1)-1;
            centreY = centreY+recording.textureRect(2)-1;
            
            [~,~,~,fileSuffix] = regexp(sta,'_channel_[0-9]+_cluster_[0-9]+.*');
            
            save([fileDir '\rf_' recording.dataFile fileSuffix{1}(1:end-4) '.mat'],'sd','rf','centreX','centreY','centreT','sds','mses');
        end
    end
end