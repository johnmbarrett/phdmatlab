function [recordingStartTime,recordingStartDate] = getMCDStartTime(file,filename)
    if isstruct(file)
        fileInfo = file;
    else
        if ischar(file)
            filename = file;
            [result,file] = ns_OpenFile(filename);
            
            if result
                error(['Could not open file: ' filename]);
            end
        elseif nargin < 2
            filename = 'unknown MCD file';
        end
        
        [result,fileInfo] = ns_GetFileInfo(file);

        if result
            error(['Could not read file info from file: ' filename]);
        end
    end

    recordingStartDate = datenum([ ...
        fileInfo.Time_Year ...
        fileInfo.Time_Month ...
        fileInfo.Time_Day ...
        fileInfo.Time_Hour ...
        fileInfo.Time_Min ...
        fileInfo.Time_Sec + fileInfo.Time_MilliSec/1000]);
    recordingStartTime = recordingStartDate*24*60*60;
end