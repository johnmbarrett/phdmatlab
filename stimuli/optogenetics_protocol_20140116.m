pws = [5 10 25 50 75 100];
nPWs = numel(pws);
reps = 10;
IPI = 2000;
outfileprefix = 'stim_file_20140116';

pwOrder = [];
for ii = 1:reps
    pwOrder = [pwOrder pws(randperm(nPWs))]; %#ok<AGROW>
end

fout = fopen(sprintf('%s_%d.txt',outfileprefix,1),'w');

for ii = 1:numel(pwOrder)
    pw = pwOrder(ii);
    
    fprintf(fout,'r 0 0 16 16 %d\r\n',pw);
    fprintf(fout,'b %d\r\n',IPI-pw);
end

fclose(fout);