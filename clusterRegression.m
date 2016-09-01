function [as,bs,xk,yk,es,k] = clusterRegression(x,y,maxK,debug)
    if numel(y) ~= numel(x);
        error('no u');
    end
    
    ess = cell(maxK,1);
    ass = cell(maxK,1);
    bss = cell(maxK,1);
    xks = cell(maxK,1);
    yks = cell(maxK,1);
        
    colours = distinguishable_colors(maxK);

    if debug
        figure;
    end

    xks{1} = {x};
    yks{1} = {y};
    coeffs = regress(y,[ones(size(x)) x]);
    ass{1} = coeffs(1);
    bss{1} = coeffs(2);
    ess{1} = sum((y - bss{1}*x - ass{1}).^2);
    
    k = 1;
    while k <= maxK
        xl = cell(0,1);
        yl = cell(0,1);
        al = [];
        bl = [];
        el = [];
        
        for ii = 1:k
            xi = xks{k}{ii};
            yi = yks{k}{ii};
            r = yi - bss{k}(ii)*xi - ass{k}(ii);
            
            [dip,p] = hartigansdipsigniftest(r,1000);
            
            if debug
                clf;
                hist(r,ceil(sqrt(numel(r))));
                
                disp(dip);
                disp(p);
            end
            
            if p < 0.05
                [as,bs,xk,yk,es] = kregress(xi,yi,2,debug,colours);
                idx = numel(xl)+(1:numel(xk));
                xl(idx) = xk;
                yl(idx) = yk;
                al(idx) = as; %#ok<*AGROW>
                bl(idx) = bs;
                el(idx) = sum(es);
            else
                xl{end+1} = xi;
                yl{end+1} = yi;
                al(end+1) = ass{k}(ii);
                bl(end+1) = bss{k}(ii);
                el(end+1) = ess{k}(ii);
            end
        end
        
        if numel(xl) == k
            break;
        end
        
        k = numel(xl);
        
        if debug
            clf;
            hold on;
            
            for ii = 1:k
                plot(xl{ii},yl{ii},'LineStyle','none','Marker','o','Color',colours(ii,:));
                plot(xl{ii},bl(ii)*xl{ii}+al(ii),'Color',1-colours(ii,:));
            end
        end
        
        ess{k} = el;
        ass{k} = al;
        bss{k} = bl;
        xks{k} = xl;
        yks{k} = yl;
    end
    
    es = ess{k}';
    as = ass{k}';
    bs = bss{k}';
    xk = xks{k};
    yk = yks{k};
end

function [as,bs,xk,yk,es] = kregress(x,y,k,debug,colours)
    n = numel(x);
    m = ceil(n/k);

    as = zeros(1,k);
    bs = zeros(1,k);

    [ys,sortIndices] = sort(y);
    xs = x(sortIndices);
    oldPartitions = zeros(size(x));

    for ii = 1:k
        idx = (ii-1)*m+1:min(ii*m,n);
        xk = xs(idx);
        yk = ys(idx);
        oldPartitions(sortIndices(idx)) = ii;
        coeffs = regress(yk,[ones(size(xk)) xk]);
        as(ii) = coeffs(1);
        bs(ii) = coeffs(2);
    end

    iterations = 0;
    while true
        es = (repmat(y,1,k) - repmat(x,1,k).*repmat(bs,n,1) - repmat(as,n,1)).^2;

        [~,newPartitions] = min(es,[],2);
        iterations = iterations + 1;

        if debug
            clf;
            hold on;

            for ii = 1:k
                plot(x(newPartitions == ii),y(newPartitions == ii),'Color',colours(ii,:),'LineStyle','none','Marker','o');
                plot(x,bs(ii)*x+as(ii),'Color',ones(1,3)-colours(ii,:));
            end

            disp(k);
            disp(iterations);
            disp(as);
            disp(bs);
            disp(sum(es))

            for ii = 1:k
                fprintf('%d\t',sum(newPartitions == ii));
            end

            fprintf('\n---\n');
        end

        if isequal(oldPartitions,newPartitions)
            break;
        end

        xk = cell(k,1);
        yk = cell(k,1);

        for ii = 1:k
            xk{ii} = x(newPartitions == ii);
            yk{ii} = y(newPartitions == ii);
            bs(ii) = corr(xk{ii},yk{ii})*std(yk{ii})/std(xk{ii});
            as(ii) = mean(yk{ii})-bs(ii)*mean(xk{ii});
        end

        oldPartitions = newPartitions;
    end
end