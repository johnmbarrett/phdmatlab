function index = getMCSSubplotIndex(chlabel)
    index = 8*(mod(chlabel,10)-1)+floor(chlabel/10);
end