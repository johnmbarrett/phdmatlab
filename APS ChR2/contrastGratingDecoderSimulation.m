function pc = contrastGratingDecoderSimulation(a,b,k,t)
    g = 5/8;
    c = 0.1:0.1:0.6;
    Imin = g*(1-c);
    Imax = g*(1+c);
    % disp(Imax);
    % disp(Imin);
    % disp((Imax-Imin)./(Imax+Imin));

    pc = zeros(6,4);

    % a = 1/(1-g);
    % b = -1*g/(1-g);

    % figure;
    % fplot(@(x) max(0,a*x+b),[0 1])

    %%

    phase = kron([1 2 3 4; 2 3 4 1; 3 4 1 2; 4 1 2 3],ones(125,1));
    si = gamrnd(k,t,500,1);

    for kk = 1:6
        xs = [Imax(kk); g; Imin(kk); g];

        for ll = 1:4
            li = si+a*xs(phase(:,ll))+b;
            logLi = log(li);
            mi = si+a*g+b;

%             ri = 0:170;

%             Li = repmat(li,1,171);
%             Mi = repmat(mi,1,171);
%             Ri = repmat(ri,500,1);

%             Pri = (Li.^Ri.*exp(-Li)+Mi.^Ri.*exp(-Mi))./(2*factorial(Ri));
%             Pri(isinf(Pri)) = 0;

%             Cri = cumsum(Pri,2);

%             Ai = log(li)-log(mi);
%             Bi = mi-li;

            pn = zeros(50,1);

            for jj = 1:50
                r1 = poissrnd(li);
                r2 = poissrnd(mi);
                d = (r1-r2);
                s = sign(-d);
                n = min(r1,r2)+1;
                m = max(r1,r2);
                LL = d.*logLi+s.*arrayfun(@(n,m) sum(log(n:m)),n,m);
                pn(jj) = sum(LL) > 0;
                continue
                
    %             tic;
                r = zeros(500,1);

                for ii = 1:500
                    u = rand(1);
                    r(ii) = ri(find(Cri(ii,:) > u,1));
            %         r(ii) = randsample(ri,1,true,Pri(ii,:));
                end

                pn(jj) = sum(Ai.*r+Bi) > 0;
    %             toc;
            end

            pc(kk,ll) = sum(pn == 1)/50;
        end

    %     figure;
    %     plot(c,pc);
    end
end