function warped = infWarp(allTimes)
    warped = allTimes;
    warped(isfinite(allTimes)) = 2*(tiedrank(allTimes(isfinite(allTimes)))-0.5)/sum(isfinite(allTimes))-1;
end