% function [ps,prk,prks] = trainDiscreteDecoder(trainIndices,s,R,vk,nk,m,N)
function pdfs = trainDiscreteDecoder(trainIndices,~,~,params)
    pdfs = params;
    
    pdfs(1).prks = cellfun(@(binr) trainDiscreteDecoderHelper(trainIndices,binr),params.binR,'UniformOutput',false);
    
%     m = params.m;
%     nk = params.nk;
%     vk = params.vk;
%     R = params.R;
%     
%     N = numel(trainIndices);
%     prks = cellfun(@(v) zeros(numel(v),m),vk,'UniformOutput',false);
%     
%     for ii = 1:m
%         sameStimulus = s(trainIndices) == ii;
%         ns = sum(sameStimulus);
% 
%         for jj = 1:numel(nk)
%             n = nk(jj);
% 
%             for kk = 1:n
%                 nrs = sum(sameStimulus & R(trainIndices,jj) == kk);
%                 prks{jj}(kk,ii) = nrs/ns;
%             end
%         end
%     end
%     
%     pdfs = struct('nk',nk,'prks',{prks},'vk',{vk});
end

function prs = trainDiscreteDecoderHelper(trainIndices,binr)
    binr = binr(trainIndices,:,:);
    
    nrs = permute(sum(binr,1),[2 3 1]);
    
    prs = nrs./repmat(sum(nrs,1),size(nrs,1),1);
end