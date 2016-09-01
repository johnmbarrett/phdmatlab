intensities = 1:2:15;
durations = 4:2:10;
reps = 10;

conditions = [kron(durations',ones(8,1)) repmat(intensities',4,1)];
conditionOrder = zeros(320,1);

for ii = 1:reps
    conditionOrder(32*(ii-1)+(1:32)) = randperm(32);
end

outfile = 'stim_file_20141107_2.txt';

fout = fopen(outfile,'w');

for ii = 1:320
    condition = conditions(conditionOrder(ii),:);
    fprintf(fout,'r 0 0 16 16 %d %d\r\n',condition(1),condition(2));
    fprintf(fout,'b %d\r\n',80-condition(1));
end

fclose(fout);
