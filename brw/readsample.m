% this sample script use the BrwExtReader.dll and allows you to load and 
% display values contained in a BRW file. Having located the library and 
% the BRW file, you can choose the channel to plot among the ones that have
% been recorded in the BRW file. Hence, the script will display, whenever
% available, both raw and spike recorded data for a maximum of 5 sec.
%
% Created on Dec 07 2011 @3Brain (info@3brain.com)

[fileName filePath] = uigetfile('*.dll','Select the BrwExtReader.dll');
libPath = strcat(filePath,fileName);
if isempty(libPath)
    return
end

[fileName filePath] = uigetfile('*.brw','Select the BRW file to load');
brwPath = strcat(filePath,fileName);
if isempty(brwPath)
    return
end

%% Open the BRW-file

% load the private assebly by the full path (comprising file name and 
% extension)
asm = NET.addAssembly(libPath);

% create an instance of the BrwFile class
brwRdr = BW.BrwRdr;

% open a BRW-file (the returned value indicates whether the opening was 
% successful)
succ = brwRdr.Open(brwPath);

if ~succ
    msgbox('Impossible to open the file. The file might be in use or the version might be not supported.')
else
    try        
        %% Get recorded channels and allow user to select one to load
        
        % gets the union of all the recorded channels regardless the 
        % recorded stream
        chs = brwRdr.GetRecChsUnion();
        
        % the number of recorded channels
        nChs = double(chs.Length);
        
        % create a list to allow user to select a channel
        chsList = cell(1, nChs);
        for i = 1 : nChs
            chsList{i} = char(chs(i).ToString());
        end
        
        [s,v] = listdlg('PromptString', 'Select a channel to load:',...
                'SelectionMode', 'single', 'ListString',chsList);
            
        % the selected channel
        ch = chs(s);        

        %% Plot data for the selected channel
                
        % sampling frequency
        sf = double(brwRdr.SamplingRate);
        
        % the number of frames to load (the min between 5 seconds and the
        % total experiment duration)
        nFrames = min(sf * 5, double(brwRdr.RecNFrames));
        
        figure
        hold on
        axis([0 nFrames/sf -4120 4120])
        title(strcat('Ch ', num2str(ch.Row), ',', num2str(ch.Col)), 'FontSize', 16)
        xlabel('Time [s]', 'FontSize', 12)
        ylabel('Values [uVolt]', 'FontSize', 12)
        
        % check whether the BRW-file contains a raw stream for the selected
        % channel
        if brwRdr.IsRawStreamFor(ch)

            % load raw data
            chYRaw = brwRdr.GetRawData(ch, 0, nFrames);

            % x values (scale in seconds)
            x = (1 : nFrames) / sf;
        
            % plot raw data
            plot(x, chYRaw, 'b');
        end

        % check whether the BRW-file contains a spike stream for the 
        % selected channel
        if brwRdr.IsSpikeStreamFor(ch)

            % load spikes
            chSpk = brwRdr.GetSpikes(ch, 0, nFrames);
            
            % the position at which spike markers will be placed on Y axis
            spkYPos = ones(1, chSpk.Length) * (-1000);
            
            % the position of spikes in the X axis (scale in seconds)
            spkXPos = double(chSpk) / sf;            

            % plot spikes
            plot(spkXPos, spkYPos, '*g');
        end
        
        % close the BRW-file
        brwRdr.Close();
        
    catch ex
        brwRdr.Close();
        throw(ex)
    end
end

