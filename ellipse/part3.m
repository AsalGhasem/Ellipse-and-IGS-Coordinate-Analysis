%% P3 - Incomplete Ellipse
close all
clear 
clc

% Final parameters (using the derived parameters from the prev part):
a  = 4.991085;
b  = 1.982314;
x0 = 1.989005;
y0 = 3.016496;

data = load('Data.mat');
XWhole = data.X_noisy;
YWhole = data.Y_noisy;

theta      = atan2((YWhole - y0)/b, (XWhole - x0)/a);  
thetaStart = deg2rad(0);              
InLane     = (theta >= thetaStart) & (theta <= thetaStart + pi/3);

X = XWhole(InLane);
Y = YWhole(InLane);

numberOfPoints = length(X);
numberOfObservations = numberOfPoints*2;
numberOfParams       = 4;

Q = eye(numberOfObservations);

syms xa ya aa bb xx0 yy0 real
F = ((xa - xx0)^2)/(aa^2) + ((ya - yy0)^2)/(bb^2) - 1;

tol = 1e-6;
maxIter = 15;
deltaNorm = inf;
iter = 0;

while deltaNorm > tol && iter < maxIter
    iter = iter + 1;
    
    J = jacobian(F, [aa, bb, xx0, yy0]);
    A = zeros(numberOfPoints, 4);
    for i = 1:numberOfPoints
        JRow = subs(J, [xa, ya, aa, bb, xx0, yy0], ...
                    [X(i), Y(i), a, b, x0, y0]);
        A(i,:) = double(JRow);
    end
    
    JB = jacobian(F, [xa, ya]);
    B = zeros(numberOfPoints, 2*numberOfPoints);
    for i = 1:numberOfPoints
        JRowB = subs(JB, [xa, ya, aa, bb, xx0, yy0], ...
                         [X(i), Y(i), a, b, x0, y0]);
        B(i, 2*i-1 : 2*i) = double(JRowB);
    end

    W = zeros(numberOfPoints,1);
    for i = 1:numberOfPoints
        Fi = subs(F, [xa, ya, aa, bb, xx0, yy0], ...
                      [X(i), Y(i), a, b, x0, y0]);
        W(i) = double(Fi);
    end
    
    P = inv(B * transpose(B));
    delta = -(inv(transpose(A) * P * A) * (transpose(A) * P * W));

    a  = a  + delta(1);
    b  = b  + delta(2);
    x0 = x0 + delta(3);
    y0 = y0 + delta(4);
    
    deltaNorm = norm(delta);

    if norm(delta) > 1000
    error("Adjustment diverging");
    end
    
    fprintf('Itertion %d: delta norm = %.6f, a=%.4f, b=%.4f, x0=%.4f, y0=%.4f\n', ...
            iter, deltaNorm, a, b, x0, y0);
end

fprintf('\nFinal parameters:\n a = %.6f\n b = %.6f\n x0 = %.6f\n y0 = %.6f\n', a, b, x0, y0);

sigmaxHat = inv(transpose(A) * P * A)
W = zeros(numberOfPoints,1);
    for i = 1:numberOfPoints
        Fi = subs(F, [xa, ya, aa, bb, xx0, yy0], ...
                      [X(i), Y(i), a, b, x0, y0]);
        W(i) = double(Fi);
    end

Ve = (A * delta) - W; % is the formula correct?
V  = Q * transpose(B) * P * Ve;

[df, testVariance, factorVariance, chi2Lower, chi2Upper] = chi2Test(V, Q, numberOfPoints, numberOfParams)

l = zeros(2*numberOfPoints,1);
for i = 1:numberOfPoints
    l(2*i-1,1) = X(i);
    l(2*i,1) = Y(i);
end
lHat = l + V;

% Residuals plot
VX = zeros(numberOfPoints,1);
VY = zeros(numberOfPoints,1);

for i = 1:numberOfPoints
    VX(i) = V(2*i-1,1);
    VY(i) = V(2*i,1);
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

% Stds plot
stds = [sqrt(sigmaxHat(1,1)); sqrt(sigmaxHat(2,2)); sqrt(sigmaxHat(3,3)); sqrt(sigmaxHat(4,4))];
params = ["a", "b", "x0", "y0"];

figure;
bar(params, stds,FaceColor = "#7E2F8E");
title('Std Of Parameters');
xlabel('Parameters');
ylabel('Std [m]');
grid on;