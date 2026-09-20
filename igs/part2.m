%% P2 - Ellipsoid
close all
clear 
clc

data = readtable("igs_coordinates.csv");

X = data.X;
Y = data.Y;
Z = data.Z;

sX = data.sX;
sY = data.sY;
sZ = data.sZ;

numberOfPoints = length(X);
numberOfObservations = numberOfPoints*3;
numberOfParams       = 5;

% Q = diag([sX.^2; sY.^2; sZ.^2]);
Q = zeros(3*numberOfPoints);
for i = 1:numberOfPoints
    Q(3*i-2,3*i-2) = sX(i)^2;
    Q(3*i-1,3*i-1) = sY(i)^2;
    Q(3*i,  3*i  ) = sZ(i)^2;
end

a  = 6378137;
b  = 6356752;
x0 = 0;
y0 = 0;
z0 = 0;

syms xa ya za aa bb xx0 yy0 zz0
F = ((xa - xx0)^2)/(aa^2) + ((ya - yy0)^2)/(aa^2) + ((za - zz0)^2)/(bb^2) - 1;

tol = 1e-4;
maxIter = 30;
deltaNorm = inf;
iter = 0;

JA = jacobian(F, [aa, bb, xx0, yy0, zz0]);
JB = jacobian(F, [xa, ya, za]);

while deltaNorm > tol && iter < maxIter
    iter = iter + 1;
    
    A = zeros(numberOfPoints, 5);
    for i = 1:numberOfPoints
        JRowA = subs(JA, [xa, ya, za, aa, bb, xx0, yy0, zz0], ...
                    [X(i), Y(i), Z(i), a, b, x0, y0, z0]);
        A(i,:) = double(JRowA);
    end
    
    B = zeros(numberOfPoints, 3*numberOfPoints);
    for i = 1:numberOfPoints
        JRowB = subs(JB, [xa, ya, za, aa, bb, xx0, yy0, zz0], ...
                         [X(i), Y(i), Z(i), a, b, x0, y0, z0]);
        B(i, 3*i-2 : 3*i) = double(JRowB);
    end

    W = zeros(numberOfPoints,1);
    for i = 1:numberOfPoints
        Fi = subs(F, [xa, ya, za, aa, bb, xx0, yy0, zz0], ...
                      [X(i), Y(i), Z(i), a, b, x0, y0, z0]);
        W(i) = double(Fi);
    end
    
    P = inv(B * Q * transpose(B));
    delta = -(inv(transpose(A) * P * A) * (transpose(A) * P * W));
   
    a  = a  + delta(1);
    b  = b  + delta(2);
    x0 = x0 + delta(3);
    y0 = y0 + delta(4);
    z0 = z0 + delta(5);
    
    deltaNorm = norm(delta);

    if norm(delta) > 1e6
    error("Adjustment diverging");
    end
    
    fprintf('Iteration %d: delta norm = %.6f, a=%.4f, b=%.4f, x0=%.4f, y0=%.4f, z0=%.4f\n', ...
            iter, deltaNorm, a, b, x0, y0, z0);
end

fprintf('\nFinal parameters:\n a = %.6f\n b = %.6f\n x0 = %.6f\n y0 = %.6f\n z0 = %.6f\n', a, b, x0, y0, z0);

sigmaxHat = inv(transpose(A) * P * A);
W = zeros(numberOfPoints,1);
    for i = 1:numberOfPoints
        Fi = subs(F, [xa, ya, za, aa, bb, xx0, yy0, zz0], ...
                      [X(i), Y(i), Z(i), a, b, x0, y0, z0]);
        W(i) = double(Fi);
    end

Ve = (A * delta) - W; % is the formula correct?
V  = Q * transpose(B) * P * Ve;

[df, testVariance, factorVariance, chi2Lower, chi2Upper] = chi2Test(V, Q, numberOfPoints, numberOfParams)

l = zeros(3*numberOfPoints,1);
for i = 1:numberOfPoints
    l(3*i-2,1) = X(i);
    l(3*i-1,1) = Y(i);
    l(3*i,1) = Z(i);
end
lHat = l + V;

% Residuals plot
VX = zeros(numberOfPoints,1);
VY = zeros(numberOfPoints,1);
VZ = zeros(numberOfPoints,1);

for i = 1:numberOfPoints
    VX(i) = V(3*i-2,1);
    VY(i) = V(3*i-1,1);
    VZ(i) = V(3*i,  1);
end

figure;
bar(1:numberOfPoints, -VX, FaceColor = "#A0F");
title('Residuals Of X');
xlabel('Point Number');
ylabel('Residual (m)');
grid on;

figure;
bar(1:numberOfPoints, -VY, FaceColor = "#0B0");
title('Residuals Of Y');
xlabel('Point Number');
ylabel('Residual (m)');
grid on;

figure;
bar(1:numberOfPoints, -VZ, FaceColor = "#00F");
title('Residuals Of Z');
xlabel('Point Number');
ylabel('Residual (m)');
grid on;

% Stds plot
stds = [sqrt(sigmaxHat(1,1)); sqrt(sigmaxHat(2,2)); sqrt(sigmaxHat(3,3)); sqrt(sigmaxHat(4,4)); sqrt(sigmaxHat(5,5))];
params = ["a", "b", "x0", "y0", "z0"];

figure;
bar(params, stds,FaceColor = "#7E2F8E");
title('Std Of Parameters');
xlabel('Parameters');
ylabel('Std [m]');
grid on;

% Observation's plot
XHat = zeros(numberOfPoints,1);
YHat = zeros(numberOfPoints,1);
ZHat = zeros(numberOfPoints,1);

for i = 1:numberOfPoints
    XHat(i) = lHat(3*i-2,1);
    YHat(i) = lHat(3*i-1,1);
    ZHat(i) = lHat(3*i,  1);
end

figure;
bar(1:numberOfPoints, XHat, FaceColor = "#A0F");
title('Observations');
xlabel('Point Number');
ylabel('X (m)');
grid on;

figure;
bar(1:numberOfPoints, YHat, FaceColor = "#0B0");
title('Observations');
xlabel('Point Number');
ylabel('Y (m)');
grid on;

figure;
bar(1:numberOfPoints, ZHat, FaceColor = "#00F");
title('Observations');
xlabel('Point Number');
ylabel('Z (m)');
grid on;

% Stds plot
stds = [sqrt(sigmaxHat(1,1)); sqrt(sigmaxHat(2,2)); sqrt(sigmaxHat(3,3)); sqrt(sigmaxHat(4,4)); sqrt(sigmaxHat(5,5))];
params = ["a", "b", "x0", "y0", "z0"];

figure;
bar(params, stds,FaceColor = "#7E2F8E");
title('Std Of Parameters');
xlabel('Parameters');
ylabel('Std [m]');
grid on;

% Ellipsoid plot
figure('Color','white')
hold on
axis equal on
grid on
view(35,25)

nLon = 120;
nLat = 60;

lon = linspace(0,2*pi,nLon);
lat = linspace(-pi/2,pi/2,nLat);
[Lon,Lat] = meshgrid(lon,lat);

Xe = x0 + a*cos(Lat).*cos(Lon);
Ye = y0 + a*cos(Lat).*sin(Lon);
Ze = z0 + b*sin(Lat);

surf(Xe,Ye,Ze,...
    'FaceColor',[0.1 0.5 0.9],...   
    'EdgeColor','none',...
    'FaceAlpha',0.9)

for L = 0:15:345
    lam = deg2rad(L);
    latLine = linspace(-pi/2,pi/2,200);

    xm = x0 + a*cos(latLine).*cos(lam);
    ym = y0 + a*cos(latLine).*sin(lam);
    zm = z0 + b*sin(latLine);

    plot3(xm,ym,zm,'k','LineWidth',0.6)
end

for B = -75:15:75
    phi = deg2rad(B);
    lonLine = linspace(0,2*pi,400);

    xp = x0 + a*cos(phi).*cos(lonLine);
    yp = y0 + a*cos(phi).*sin(lonLine);
    zp = z0 + b*sin(phi)*ones(size(lonLine));

    plot3(xp,yp,zp,'k','LineWidth',0.6)
end

scatter3(X,Y,Z,8,'r','filled')
camlight headlight
lighting gouraud
material shiny

title('Adjusted Earth Ellipsoid','FontSize',14)