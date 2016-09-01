function nsObjInfo = mungeFINDToNSObj(findInfo,type)
    nsObjInfo = struct(zeros(numel(findInfo),0));
    
    for ii = 1:numel(findInfo)
        switch(type)
            case 'FileInfo'
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
            case 'SegmentInfo'
                nsObjInfo(ii).dSampleRate           = findInfo(ii).SampleRate;
                nsObjInfo(ii).szUnits               = 'V'; % TODO : or is it?
            case 'SegmentSourceInfo'
                nsObjInfo(ii).dMinVal               = findInfo(ii).MinVal;
                nsObjInfo(ii).dMaxVal               = findInfo(ii).MaxVal;
                nsObjInfo(ii).dResolution           = findInfo(ii).Resolution;
                nsObjInfo(ii).dSubSampleShift       = findInfo(ii).SubSampleShift;
                nsObjInfo(ii).dLocationX            = findInfo(ii).LocationX;
                nsObjInfo(ii).dLocationY            = findInfo(ii).LocationY;
                nsObjInfo(ii).dLocationZ            = findInfo(ii).LocationZ;
                nsObjInfo(ii).dLocationUser         = findInfo(ii).LocationUser;
                nsObjInfo(ii).dHighFreqCorner       = findInfo(ii).HighFreqCorner;
                nsObjInfo(ii).dwHighFreqOrder       = findInfo(ii).HighFreqOrder;
                nsObjInfo(ii).szHighFilterType      = findInfo(ii).HighFilterType;
                nsObjInfo(ii).dLowFreqCorner        = findInfo(ii).LowFreqCorner;
                nsObjInfo(ii).dwLowFreqOrder        = findInfo(ii).LowFreqOrder;
                nsObjInfo(ii).szLowFilterType       = findInfo(ii).LowFilterType;
                nsObjInfo(ii).szProbeInfo           = findInfo(ii).ProbeInfo;
        end
    end
end