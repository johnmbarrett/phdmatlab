A1 = [3 2 3 3 2 3 3 3];
B1 = [2 1 2 2 1 2 2 2];
x1 = [0 0 0 0 0 0 0 0];
y1 = [0 0 0 0 0 0 0 0];
p1 = [0 0 0 0 0 0 0 0];

A2 = [2 1 2 3 2 2 1 3];
B2 = [1 0.75 1 3 1 1 2 2];
x2 = [-2 -2.5 -0.75 0 0 -1.0245209260022 0 0];
y2 = [-1 1.5 0.25 1 2 0.25 0 0];
p2 = [pi/4 pi/4 pi/4 0 0 pi/4 0 0];

%%

intersectionArea = [3.82254574 0.00000000 6.28318531 17.60218852 0.00000000 6.28318531 6.285318531 pi*A1(8)*B1(8)];
unionArea = pi*(A1.*B1+A2.*B2)-intersectionArea;
areaRatio = intersectionArea./unionArea;

%%

[xA,yA] = ellipseFn(x1(1),y1(1),A1(1),B1(1),p1(1));
[xB,yB] = ellipseFn(x2(1),y2(1),A2(1),B2(1),p2(1));

figure;
hold on;
ezplot(xA,yA,[-pi pi]);
ezplot(xB,yB,[-pi pi]);

%%

alphax = atan2([B1;B2].*sin([p1;p2]),[A1;A2].*cos([p1;p2]));
alphay = atan2([B1;B2].*cos([p1;p2]),[A1;A2].*sin([p1;p2]));

thetax = [pi-alphax -alphax];
thetay = [pi+alphay  alphay];

xlims = [x1 x1; x2 x2]+[A1 A1; A2 A2].*cos(thetax).*cos([p1 p1; p2 p2])-[B1 B1; B2 B2].*sin(thetax).*sin([p1 p1; p2 p2]);
ylims = [y1 y1; y2 y2]+[A1 A1; A2 A2].*cos(thetay).*sin([p1 p1; p2 p2])+[B1 B1; B2 B2].*sin(thetay).*cos([p1 p1; p2 p2]);

xmins = floor(min(xlims(:,1:8)));
ymins = floor(min(ylims(:,1:8)));

xmaxs = ceil(max(xlims(:,9:16)));
ymaxs = ceil(max(ylims(:,9:16)));

%%

e = 0:0.25:3.5;
n = zeros(15,8);
t = zeros(15,8);
I = zeros(15,8);
U = zeros(15,8);
R = zeros(15,8);

for ii = 1:15
    delta = 10^(-e(ii));
    
    for jj = 1:8
        xmin = xmins(jj);
        xmax = xmaxs(jj);
        ymin = ymins(jj);
        ymax = ymaxs(jj);
        
        tic;
        [Y,X] = ndgrid(ymin:delta:ymax,xmin:delta:xmax);
        n(ii,jj) = numel(Y);
        
        E1 = ((X-x1(jj))*cos(p1(jj))+(Y-y1(jj))*sin(p1(jj))).^2/A1(jj)^2+((Y-y1(jj))*cos(p1(jj))-(X-x1(jj))*sin(p1(jj))).^2/B1(jj)^2 <= 1;
        E2 = ((X-x2(jj))*cos(p2(jj))+(Y-y2(jj))*sin(p2(jj))).^2/A2(jj)^2+((Y-y2(jj))*cos(p2(jj))-(X-x2(jj))*sin(p2(jj))).^2/B2(jj)^2 <= 1;
        R(ii,jj) = sum(sum(E1 & E2))/sum(sum(E1 | E2));
        t(ii,jj) = toc;
        fprintf('Calc ratio for pair %d at res %d in %f seconds\n',jj,ii,t(ii,jj));
        I(ii,jj) = sum(sum(E1 & E2))*delta^2;
        U(ii,jj) = sum(sum(E1 | E2))*delta^2;
    end
end

%%

dI = abs(I-repmat(intersectionArea,15,1));
dU = abs(U-repmat(unionArea,15,1));
dR = abs(R-repmat(areaRatio,15,1));

pI = 100*dI./repmat(intersectionArea,15,1);
pU = 100*dU./repmat(unionArea,15,1);
pR = 100*dR./repmat(areaRatio,15,1);

%%

close all;

%%

x = reshape(n(:,[1 3 4 6 7 8]),15*6,1); 
y = reshape(pI(:,[1 3 4 6 7 8]),15*6,1);

figure;
scatter(n(:),pI(:));
set(gca,'XScale','log');
xlabel('n');
ylabel('% Error');

%%

figure;
scatter(n(:),t(:));
set(gca,'XScale','log');
xlabel('n');
ylabel('Time (seconds)');

%%

fun = @(b,x) b(1).*exp(-x.*b(2));
beta = nlinfit(x,y,fun,[1 1]);
hold on;
fplot(@(x) fun(beta,x),xlim,'Color','r');

%%

% n = ((xmax-xmin)*10.^e+1).*((ymax-ymin)*10.^e+1);
figure;
plot(n,[dI dU dR]);
set(gca,'XScale','log','XTick',n);

figure;
plot(n,[pI pU pR]);
set(gca,'XScale','log','XTick',n);

figure;
plot(n,t);
set(gca,'XScale','log','XTick',n);