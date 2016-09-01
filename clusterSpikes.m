function clusterSpikes(recording,varargin)

	% CLUSTERMYSPIKES   Cluster spike that have been detected by GetMeMySpikes
	%
	%	SYNTAX
	%		 = CLUSTERMYSPIKES(I_RECORDINGS,VARARGIN)
	%
	%   ARGUMENTS
	%		I_RECORDINGS =  
	%		VARARGIN =  
	%
	%   DESCRIPTION
	%       Long description
	%

	%
	%	Created by Jonas Zimmermann on 2009-06-11.
	%	Copyright (C)  Newcastle University. All rights reserved.
	%
    
    currentDir = pwd;
    restoreDir = onCleanup(@() cd(currentDir));

    options = getopt('fetchremote=''no'' channel chlabel pushremote=''no'' overwrite=''no''',varargin);
    if isempty(options.channel)
		options.channel = [];
    end
    if isempty(options.chlabel)
		options.chlabel = [];
    end

    parameters = struct();
	parameters.w_pre=30;                       %number of pre-event data points stored
	parameters.w_post=34;                      %number of post-event data points stored
	parameters.detection = 'neg';              %type of threshold
%	parameters.stdmin = 5.400;                  %minimum threshold
    	parameters.stdmin = 5.000;                  %minimum threshold %HDK
	parameters.stdmax = 65;                    %maximum threshold
	parameters.interpolation = 'y';            %interpolation for alignment
	parameters.int_factor = 4;                 %interpolation factor
	parameters.detect_fmin = 300;              %high pass filter for detection (default 300)
	parameters.detect_fmax = 3000;             %low pass filter for detection (default 3000)
	parameters.sort_fmin = 300;                %high pass filter for sorting (default 300)
	parameters.sort_fmax = 3000;               %low pass filter for sorting (default 3000)

	%parameters.max_spk = 10000;                 %C.E.
    parameters.max_spk = 8000;                % max. # of spikes before starting templ. match.
	parameters.template_type = 'center';       % nn, center, ml, mahal
	parameters.template_sdnum = 3;             % max radius of cluster in std devs.
	%parameters.permut = 'y';                   % for selection of random 'par.max_spk' spikes before starting templ. match. 
	 parameters.permut = 'n';                 % for selection of the first 'par.max_spk' spikes before starting templ. match.

	parameters.features = 'wav';               %choice of spike features
	parameters.inputs = 10;                    %number of inputs to the clustering
	parameters.scales = 4;                     %scales for wavelet decomposition
	if strcmp(parameters.features,'pca');      %number of inputs to the clustering for pca
		parameters.inputs=3; 
	end

	parameters.mintemp = 0;                    %minimum temperature
	parameters.maxtemp = 0.301;                %maximum temperature
	parameters.tempstep = 0.01;                %temperature step
	parameters.num_temp = floor(...
	(parameters.maxtemp - ...
	parameters.mintemp)/parameters.tempstep); %total number of temperatures 
	parameters.stab = 0.8;                     %stability condition for selecting the temperature
	parameters.SWCycles = 150;                 %number of montecarlo iterations
	parameters.KNearNeighb = 11;               %number of nearest neighbors
	parameters.randomseed = 0;                 % if 0, random seed is taken as the clock value
	%parameters.randomseed = 147;              % If not 0, random seed   
	parameters.min_clus = 10;                   % minimun size of a cluster
	parameters.max_clus = 12;                    % maximum number of clusters allowed

	parameters.min_clus_abs = 10;              %minimum cluster size (absolute value)
	parameters.min_clus_rel = 0.005;           %minimum cluster size (relative to the total nr. of spikes)
	%parameters.temp_plot = 'lin';               %temperature plot in linear scale
	parameters.temp_plot = 'log';              %temperature plot in log scale
	parameters.force_auto = 'y';               %automatically force membership if temp>3.
	parameters.max_spikes = 5000;              %maximum number of spikes to plot.

	%parameters.sr = 24000;                     %sampling frequency, in Hz.


	%%%%%%files = textread('Files.txt','%s');
	%%%%%%
	%%%%%%for k=1:length(files)
	%%%%%%    tic
	%%%%%%    dataFilePrefix = files(k);

	fig = figure;
    closeFigure = onCleanup(@() close(fig));
    
    dirName = getAnalysisOutputDir(recording);

    parameters.sr = recording.sampleRate;                     %sampling frequency, in Hz.filename = recording.dataFile;
    filename = recording.dataFile;

    if isfield(recording,'spont') 
        spont = recording.spont;
    else
        spont = false;
    end

    channelInfo = getMetaInfo(recording.index,'spont',spont);

    channels = parseChannelOptions(recording.index, 'channel',options.channel,'chlabel',options.chlabel,'spont',spont);

    for ii = 1:numel(channels);
        tic
        channel = channels(ii);
        channelLabel = channelInfo(channel).label;

        if any(spont)
            dataFile = [filename '_channel_' channelLabel];
        else
            dataFile = [filename '_cleaned_channel_' channelLabel];
        end
        
        dataFilePrefix = [dirName '\' dataFile];

        set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]) 

        try
            % load spikes
            jonasSpikeFile = [dataFile '_spikes.mat'];
            mcdSpikeFile = [dataFile '_MCD_spikes.mat'];
            trimmedSpikeFile = [dataFile '_MCD_trimmed_spikes.mat'];

            doTrim = false;
            if exist([dirName '\' jonasSpikeFile],'file');
                spikeFile = jonasSpikeFile;
            elseif exist([dirName '\' trimmedSpikeFile],'file');
                spikeFile = trimmedSpikeFile;
            elseif exist([dirName '\' mcdSpikeFile],'file');
                spikeFile = mcdSpikeFile;
                doTrim = true;
            else
                warning('No spike file found for channel %s (%d/%d) in file %d (%s)\n',channelLabel,ii,numel(channels),recording.index,recording.dataFile); %#ok<WNTAG>
                continue;
            end

            load([dirName '\' spikeFile]);
            
            if doTrim
                spikes = spikes(:,10:57);
            end

            nSpikes = size(spikes,1);
            maxNSpikes = min(parameters.max_spk,nSpikes);
            parameters.min_clus = max(parameters.min_clus_abs,parameters.min_clus_rel*maxNSpikes);

            % CALCULATES INPUTS TO THE CLUSTERING ALGORITHM. 
            inspk = wave_features(spikes,struct('par',parameters));              %takes wavelet coefficients.

            % SELECTION OF SPIKES FOR SPC 
            if parameters.permut == 'n'
                ipermut = 1:size(inspk,1);
            else
                ipermut = randperm(size(inspk,1));
            end

            % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
            if nSpikes > parameters.max_spk;
                ipermut = ipermut(1:maxNSpikes);
            end

            inputSpikes = inspk(ipermut,:); %#ok<NASGU>

            %INTERACTION WITH SPC
            cd(dirName); % Cluster.exe doesn't seem to like accessing files in other directories
            
            parameters.fname_in=['tmp_data' '_' channelLabel];
            parameters.fname = ['data_' channelLabel];   %filename for interaction with SPC

            save(parameters.fname_in,'inputSpikes','-ascii');
            [clu, tree] = run_cluster(struct('par',parameters));
            cd ..;
            
            [temp] = find_temp(tree,struct('par',parameters));

            %DEFINE CLUSTERS
            classes = zeros(nSpikes,1);

            for jj = 1:5
                classes(ipermut(clu(temp,3:end)==jj-1)) = jj;
            end

            % IF TEMPLATE MATCHING WAS DONE, THEN FORCE
            if nSpikes > parameters.max_spk || parameters.force_auto == 'y'
                validSpikes = spikes(classes~=0,:);
                noiseSpikes = spikes(classes==0,:);
                classes(classes==0) = force_membership_wc(validSpikes, classes(classes~=0,:), noiseSpikes, struct('par',parameters));
            end    
            
            % TODO : necessary?
            class = cell(1,6);
            for jj = 0:5
                class{jj+1} = find(classes == jj);
            end

            %PLOTS
            clf;
            ylimit = [];
            subplot(3,5,11);
            temperature=parameters.mintemp+temp*parameters.tempstep;
            
            plotfun = @plot;
            
            switch parameters.temp_plot
                case 'lin'
                    plotfun = @plot;
                case 'log'
                    plotfun = @semilogy;
            end
            
            hold on;
            
            minTemp = parameters.mintemp;
            nTemp = parameters.num_temp;
            tempStep = parameters.tempstep;
            tempDiff = parameters.maxtemp - tempStep;
            
            plotfun([minTemp tempDiff],repmat([parameters.min_clus],1,2),'k:');
            plotfun(minTemp+(1:nTemp)*tempStep, tree(1:nTemp,5:size(tree,2)));
            plotfun([temperature temperature],[1 tree(1,5)],'k:');
            
            subplot(3,5,6);
            hold on;
            
            cluster_class=zeros(nSpikes,2);
            cluster_class(:,2)= index';
            nSpikesPerCluster = length(class{1});
            colours = 'kbrgcmy';
            
            for jj = 1:6
                cluster = class{jj};
                
                if length(cluster) < parameters.min_clus;
                    continue;
                end
                
                if jj > 1
                    nSpikesPerCluster = [nSpikesPerCluster length(cluster)]; %#ok<AGROW>
                    cluster_class(cluster,1) = jj-1;
                end
                
                subplot(3,5,6);
                maxNSpikes = min(length(cluster),parameters.max_spikes);
                plot(spikes(cluster(1:maxNSpikes),:)',colours(jj)); 
                xlim([1 size(spikes,2)]);
                
                if jj < 5
                    subplot(3,5,5+jj+(jj==1)*4);
                    hold on
                    plot(spikes(cluster(1:maxNSpikes),:)',colours(jj));  
                    plot(mean(spikes(cluster,:),1),'Color',char((jj==1)*'c'+(jj~=1)*'k'),'linewidth',2)
                    xlim([1 size(spikes,2)]); 
                    title(['Cluster ' num2str(jj-1)],'Fontweight','bold')
                    
                    if jj > 1
                        ylimit = [ylimit;ylim]; %#ok<AGROW>
                    end
                    
                    subplot(3,5,10+jj+(jj==1)*4);
                    isis = diff(index(cluster));
                    [frequency,bins] = hist(isis,0:1:100);
                    bar(bins(1:end-1),frequency(1:end-1));
                    xlim([0 100])
                    xlabel([num2str(sum(frequency(1:3))) ' in < 3ms'])
                    title([num2str(length(cluster)) ' spikes']);
                end
            end

            % Rescale spike's axis 
            if ~isempty(ylimit)
                ymin = min(ylimit(:,1));
                ymax = max(ylimit(:,2));
            else
                ymin = -200;
                ymax = 200;
            end
            
            for jj = 1:4
                if length(class{jj}) > parameters.min_clus
                    subplot(3,5,5+jj+(jj==1)*4);
                    ylim([ymin ymax]);
                end
            end

            subplot(3,1,1)
            box off;
            hold on;
            
            %% these lines are for plotting continuous data; CPG 3/4/2008
            concatFile = [dataFilePrefix '_concat.mat'];
            rawFile = [dataFilePrefix '.mat'];
            if exist(concatFile,'file')
                load(concatFile);
            elseif exist(rawFile,'file')
                load(rawFile);
                data = channelData{1}; %#ok<USENS>
                clear channelData;
            else
                if ii == 1
                    safeLoadMCDLibrary;
                    [~,file] = ns_OpenFile([recording.dataFile '.mcd']);
                    closeDataFile = onCleanup(@() ns_CloseFile(file));
                end

                [~,~,data] = ns_GetAnalogData(file,channelInfo(channel).numInFile,1,60*parameters.sr);
            end

            if length(data)>60*parameters.sr
                data = data(1:60*parameters.sr)';       %will plot just 60 sec.
            end

            %Filters and gets threshold
            [b,a]=ellip(2,0.1,40,[parameters.sort_fmin parameters.sort_fmax]*2/(parameters.sr));
            filteredData=filtfilt(b,a,data);
            minThreshold = parameters.stdmin * median(abs(filteredData))/0.6745;
            maxThreshold = parameters.stdmax * median(abs(filteredData))/0.6745;
            
            t = (1:length(filteredData))/parameters.sr;
            plot(t,filteredData)
            
            if ~strcmp(parameters.detection,'neg')
                line([0 t(end)],[minThreshold minThreshold],'color','r')
            end
            
            if ~strcmp(parameters.detection,'pos')
                line([0 t(end)],[-minThreshold -minThreshold],'color','r')
            end
            
            ylim([-maxThreshold maxThreshold])
            
            % end of continuous data plotting; CPG 3/4/2008
            title(dataFile,'Interpreter','none','Fontsize',14)
            features = parameters.features;

            set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
            set(gcf,'paperposition',[.25 .25 10.5 7.8])
            print('-dpng','-r500',[dirName '\fig_' char(dataFile) '_cluster']);

            %SAVE FILES
            par = parameters; %#ok<NASGU> preserving the name par in the output for compatability but I refuse to use it anywhere else in the file
            outfile=[dirName '\times_' spikeFile];
            if parameters.permut == 'n'
                save(outfile, 'cluster_class', 'par', 'spikes', 'inspk');
            else
                save(outfile, 'cluster_class', 'par', 'spikes', 'inspk', 'ipermut');
            end

            nClusters=length(nSpikesPerCluster)-1;
            
            clusterResultsFile=[dataFilePrefix 'cluster_results.txt'];
            fout=fopen(clusterResultsFile,'at+');
            closeResultsFile = onCleanup(@() fclose(fout));
            
            fprintf(fout,'%s\t %s\t %g\t %d %g\t', char(dataFilePrefix), features, temperature, nClusters, parameters.stdmin);
            for jj = 1:nClusters+1
                fprintf(fout,'%d\t',nSpikesPerCluster(jj));
            end
            fprintf(fout,'%d\n',nSpikesPerCluster(end));

            fprintf('Done clustering spikes for %i (%s), channel %s (%i/%i) after %gs\n', recording.index, filename, channelLabel, ii, length(channels), toc);
        catch channelError
            logMatlabError(channelError);
            cd(currentDir);
            fprintf('Clustering spikes failed for file %i (%s), channel %s (%i/%i) after %gs\n', recording.index, filename, channelLabel, ii, length(channels), toc);
        end
    end
end