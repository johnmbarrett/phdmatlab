function pdfs = trainKDEDecoder(trainIndices,s,R,params)
%     nb = NaiveBayes.fit(R(trainIndices,:),s(trainIndices),'Distribution','kernel','KSSupport',[lb ub]);
    n = 2^12;
    k = size(R,2);
    
    if n == k % this would break makeSubDecoder but is also extremely unlikely to happen
        n = 2^13;
    end
    
    m = max(s);
    
    prks = zeros(n,k,m);
    bw = params.bw;
    lb = params.lb;
    ub = params.ub;
    x = params.x;
    
    for ii = 1:k
        for jj = 1:m
            idx = s(trainIndices) == jj;
            prks(:,ii,jj) = kernelPDFEstimate(R(trainIndices(idx),ii),x,bw(ii,jj),lb,ub);
        end
    end
    
    assert(all(isfinite(x)));
    
%     f1 = figure;
%     hold on;
%     
%     [h,b] = hist(R(trainIndices,:));
%     bar(b,h./repmat(sum(h),10,1));
%     plot(x,prk);
%     
%     f2 = figure;
%     
%     for ii = 1:m
%         subplot(2,4,ii);
%         hold on;
%         idx = s(trainIndices) == ii;
%         [hs,bs] = hist(R(trainIndices(idx),:));
%         bar(bs,hs./repmat(sum(hs),10,1));
%         plot(x,prks(:,:,ii));
%     end
%     
%     close(f1);
%     close(f2);
    
    pdfs = struct('prks',prks,'x',params.x);
end