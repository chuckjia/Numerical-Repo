clear; clc; close all

Nx = 200;

% Mesh and physical parameters
Np = Nx;
xf = 75000;

% Run the projection method
[u_after, u_before, der_after, der_before] = main(Nx);

% Plot the results
xVec = linspace(0, xf, Nx);
plot(xVec, u_after(:, Np), '.');
hold on
plot(xVec, u_before(:, Np), '.');
legend('After', 'Before')
plot(4.1438e4, 10, 'rd')
hold off

xVec = xVec(1:end-1);
figure
plot(xVec, der_after, '.')
hold on
% plot(xVec, der_before, '.')
hold off




