function [ifr,indices] = IFR(times,bin,usePauluisBakerMethod,a,alpha)
    if nargin < 2
        bin = 1;
    end
    
    indices = ceil(times/bin);
    isi = diff(times);
    
    ifr = zeros(max(indices),1);
    
    for ii = 1:numel(indices)-1
        ifr(indices(ii):indices(ii+1)-1) = 1./isi(ii);
    end
    
    if nargin < 3 || ~all(logical(usePauluisBakerMethod))
        return;
    end
    
    if nargin < 4
        phat = gamfit(isi);
        a = phat(1);
    end
    
    if nargin < 5
        alpha = 0.05;
    end
    
    increases = zeros(size(isi));
    decreases = zeros(size(isi));
    
    for ii = 1:numel(isi)
        if ii > 2
            p1 = gammaPVal(isi(ii),isi(ii-1),a);
            p2 = gammaPVal(isi(ii),isi(ii-2),a);
            decreases(ii) = (p1 < alpha) && (p2 < alpha);
        end
        
        if ii < numel(isi)-1
            p1 = gammaPVal(isi(ii),isi(ii+1),a);
            p2 = gammaPVal(isi(ii),isi(ii+2),a);
            increases(ii) = (p1 < alpha) && (p2 < alpha);
        end
    end
    
    changess = {find(increases) find(decreases)};
    
    for jj = [1 -1]
        changes = changess{(3-jj)/2};
        
        for ii = 1:numel(changes)
            idx = changes(ii);
            In = isi(idx);
            Ishort = isi(idx+jj);
            Thigh = 0.5*Ishort;
            Ihigh = ceil(Thigh/bin);
            Fhigh = 1/Ishort;
            Tlow = In - Thigh;
            Ilow = ceil(Ilow/bin);
            Flow = 1/(2*Tlow);

            ifr(indices(idx)+jj*(-Ilow:Ihigh)) = Flow;
        end
    end
end

function p = gammaPVal(isi,mu,a)
    p = 1-gamcdf(isi,a,mu/a);
end