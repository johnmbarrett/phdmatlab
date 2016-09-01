function nsObjInfo = mungeFINDFileInfoToNSObjFileInfo(findInfo)
    nsObjInfo = struct(zeros(numel(findInfo),0));
    
    for ii = 1:numel(findInfo)
        nsObjInfo(ii).szFileType            = findInfo(ii).FileType;
        nsObjInfo(ii).dTimeStampResolution  = findInfo(ii).TimeStampResolution;
        nsObjInfo(ii).dTimeSpan             = findInfo(ii).TimeSpan;
        nsObjInfo(ii).szAppName             = findInfo(ii).AppName;
        nsObjInfo(ii).dwTime_Year           = findInfo(ii).Time_Year;
        nsObjInfo(ii).dwTime_Month          = findInfo(ii).Time_Month;
        nsObjInfo(ii).dwTime_Day            = findInfo(ii).Time_Day;
        [~,recordingStartDate]              = getMCDStartTime(findInfo);
        nsObjInfo(ii).dwTime_DayOfWeek      = weekday(recordingStartDate)-1;
        nsObjInfo(ii).dwTime_Hour           = findInfo(ii).Time_Hour;
        nsObjInfo(ii).dwTime_Min            = findInfo(ii).Time_Min;
        nsObjInfo(ii).dwTime_Sec            = findInfo(ii).Time_Sec;
        nsObjInfo(ii).dwTime_MilliSec       = findInfo(ii).Time_MilliSec;
        nsObjInfo(ii).szFileComment         = findInfo(ii).FileComment;
    end
end