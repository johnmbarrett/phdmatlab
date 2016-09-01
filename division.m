I = imread('path.gif');
%%
Z = double(I < 255);
ind = find(Z);
[x,y] = ind2sub(size(Z),ind);
[xc,yc,r] = circfit(x,y);
% umppx = 0.645089285714286;
% L = 3159.847/umppx; % substitute the path length in um where it says <path length>
% a = L/r;
u = x-xc;
v = y-yc;
[theta,rho] = cart2pol(u,v);
a = max(diff(sort(theta)));
b = a/7;
%%
figure;
hold on;
minX = Inf;
maxX = -Inf;
for ii = 1:7
    theta = b*(ii-1);
    x0 = xc;
    x1 = 2*r*cos(theta)+xc;
    maxX = max(maxX,max(x0,x1));
    minX = min(minX,min(x0,x1));
    fun = @(x) tan(theta)*(x-xc)+yc;
    fplot(fun,[x0 x1]);
end
xlim([minX maxX]);
%%
[Y,X] = ndgrid(1:size(Z,1),1:size(Z,2));
Y = Y-yc;
X = X-xc;
T = atan2(Y,X);
R = sqrt(X.^2 + Y.^2);
D = zeros(size(Z));

for ii = 1:7
    theta = b*(ii-1);
    tol = min(pi/100,pi./R);
    D(abs(T-theta) < tol) = 255;
end

imwrite(flipud(255-D),'division.gif');
