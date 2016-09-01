function err = safeLoadMCDLibrary
    mcdLibrarySuffix = '\phd\matlab\nsMCDLibrary\Matlab\Matlab-Import-Filter\Matlab_Interface\nsMCDLibrary64.dll';
    serverMCDLibrary = 'V:\retina\Shared\Matlab\nsMCDLibrary64.dll';
    roamingMCDLibrary = ['H:' mcdLibrarySuffix];
    localMCDLibrary = ['D:\backup' mcdLibrarySuffix];

    if exist(localMCDLibrary,'file')
        err = ns_SetLibrary(localMCDLibrary);
    elseif exist(roamingMCDLibrary,'file')
        err = ns_SetLibrary(roamingMCDLibrary);
    elseif exist(serverMCDLibrary,'file')
        err = ns_SetLibrary(serverMCDLibrary);
    else
        err = -1; % TODO : undocumented?
    end
end