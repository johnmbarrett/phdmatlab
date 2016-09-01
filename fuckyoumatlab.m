for ii = 1:numel(spikeFiles)
load(spikeFiles{ii});
clusters = unique(cluster_class(:,1));
for jj = 1:numel(clusters)
ifr = 1./diff(cluster_class(cluster_class(:,1) == clusters(jj),2));
hist(ifr,100);
% input('...')
end
end