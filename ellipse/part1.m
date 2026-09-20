%% P1 - Normal Ellipse
close all
clear 
clc

data = load('Data.mat');
X = data.X_noisy;
Y = data.Y_noisy;

numberOfPoints = length(X);
numberOfObservations = numberOfPoints*2;
numberOfParams       = 3;

a  = 5;
b  = 2;
x0 = 2;
y0 = 2.75;

Q = eye(numberOfObservations);

syms xa ya aa bb xx0 yy0 real
F = ((xa - xx0)^2)/(aa^2) + ((ya - yy0)^2)/(bb^2) - 1;

tol = 1e-6;
maxIter = 100;
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

% Observation's plot
XHat = zeros(numberOfPoints,1);
YHat = zeros(numberOfPoints,1);

for i = 1:numberOfPoints
    XHat(i) = lHat(2*i-1,1);
    YHat(i) = lHat(2*i,1);
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

% Fitted Ellipse
figure
scatter(X,Y,50,"MarkerEdgeColor","b", ...
    "MarkerFaceColor",[0 0.7 0.7], 'DisplayName','Noisy points')
hold on
axis equal
grid on

theta = linspace(0, 2*pi, 100);
x_ellipse = x0 + a * cos(theta);
y_ellipse = y0 + b * sin(theta);
plot(x_ellipse, y_ellipse, 'black-', 'LineWidth',2, 'DisplayName','Fitted ellipse');

legend
xlabel('X[m]')
ylabel('Y[m]')
xlim([-4,8])
title('Noisy points and Fitted Ellipse using Least Squares')