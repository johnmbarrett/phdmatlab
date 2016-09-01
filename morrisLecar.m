C = 1e-6;
gL = 3.5e-5;
gCa = 5.25e-4;
gK = 1e-3;
VL = -6.5e-2;
VCa = 4e-2;
VK = -0.1;
V1 = -1.2e-3;
V2 = 2.05e-2;
V3 = 2e-3;
V4 = 1.5e-2;
phi = (1+sqrt(5))/2;
Mss = @(V) (1+tanh((V-V1)/V2))/2;
Nss = @(V) (1+tanh((V-V3)/V4))/2;
tauN = @(V) 1/(phi*cosh((V-V3)/(2*V4)));

timestep = 0.001;
maxT = 1;
t = 0;

T = t:timestep:maxT;
N = zeros(ceil(maxT/timestep)+1,1);
V = zeros(ceil(maxT/timestep)+1,1);
V(1) = VL;

for ii = 1:numel(T)-1
    t = T(ii);
    
    if t >= 0.1 && t < 0.2
        I = 1e-6;
    else
        I = 0;
    end
    
    v = V(ii);
    n = N(ii);
    dN = @(~,N) (N-Nss(v))/tauN(v);
    tic;
    [~,Y] = ode45(dN,[t T(ii+1)],N(ii));
    toc;
    N(ii+1) = Y(end);
    dV = @(~,V) (I-gL*(V-VL)-gCa*Mss(V-VCa)-gK*n*(V-VK))/C;
    tic;
    [~,Y] = ode45(dV,[t T(ii+1)],V(ii));
    toc;
    V(ii+1) = Y(end);
end