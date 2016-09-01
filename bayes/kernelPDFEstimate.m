function pdf = kernelPDFEstimate(sample,x,bw,lb,ub)
    if all(sample < lb | sample > ub)
        pdf = zeros(size(x));
        return;
    end
    
    % this is ripped wholesale from kde.m with some renaming of
    % variables
    % histogram
    h = histc(sample,x);
    h = h(:)/sum(h);

    % dct1d
    n = numel(x);
    w = [1;2*(exp(-1i*(1:n-1)*pi/(2*n))).'];
    a = [ h(1:2:end,:); h(end:-2:2,:) ];
    a = real(w.* fft(a));

    % idct1d
    t_star = (bw/(ub-lb))^2;
    b = a.*exp(-(0:n-1)'.^2*pi^2*t_star/2);
%             n = size(b,1);
    w = n*exp(1i*(0:n-1)*pi/(2*n)).';
    b = real(ifft(w.*b));

    pdf = zeros(n,1);
    pdf(1:2:n) = b(1:n/2);
    pdf(2:2:n) = b(n:-1:n/2+1);
            
    if any(pdf < 0)
        pdf = pdf-min(pdf);
    end
    
    pdf(pdf == 0) = eps;
    
    assert(all(pdf > 0));
    
    pdf = pdf/sum(pdf);
end