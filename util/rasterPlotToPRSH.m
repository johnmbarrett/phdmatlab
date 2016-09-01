% function [rasters,edgess,histograms,repeats] = rasterPlotToPRSH(spikeTimes,stimulusTimes,ts,bw,conditions,groups,factors,levels,doFigs,rasters,edgess,histograms,repeats)
%     [figs,lines,edges,hists,trials] = rasterPlot(spikeTimes,stimulusTimes,conditions,ts,[],true,bw,groups,factors);
% 
%     if isempty(rasters)
%         rasters = cell([0 levels 2]);
%         histograms = cell([0 levels 2]);
%         edgess = cell([0 levels 2]);
%         repeats = zeros([0 levels 2]);
%     end
% 
%     rasters = [rasters; cat(ndims(rasters),reshape(lines,[1 size(lines)]),cell([1 size(lines)]))];
%     edgess = [edgess; cat(ndims(edgess),reshape(repmat({edges},levels),[1 levels]),cell([1 levels]))];
% 
%     dimensionDistributions = cell(1,numel(factors));
% 
%     for ii = 1:numel(factors)
%         dimensionDistributions{ii} = ones(levels(ii),1);
%     end
% 
%     histSize = size(hists);
%     hists = [hists; zeros([1 histSize(2:end)])]; % cpsrh makes histograms with the same number of bins as edges, whereas rasterPlot has one less bin.  Hence pad out the extra bin with zeros (as rasterPlot will never put anything in the last bin)
%     histograms = [histograms; cat(ndims(histograms),mat2cell(hists,numel(edges),dimensionDistributions{:}),cell([1 levels]))];
% 
%     repeats = [repeats; cat(ndims(repeats),reshape(trials,[1 size(trials)]),zeros([1 size(trials)]))];
% 
%     if nargin < 8 || doFigs
%         return;
%     end
% 
%     figPrefix = sprintf('%s\\%s_raster_%s_channel_%s_cluster_%d',stimDir,prefix,stimFile,channel,cluster);
% 
%     for ii = 1:numel(figs)
%         set(figs(ii),'Position',[0 0 1600 900]); %,'Visible','on');
% 
%         figFile = figPrefix;
%         for jj = 2:nFields
%             figFile = sprintf('%s_%s_%d',figPrefix,factors{jj},uniqueGroups(ii,jj-1));
%         end
% 
%         saveas(figs(ii),figFile,'fig');
%         saveas(figs(ii),figFile,'png');
%         close(figs(ii));
%     end
% end