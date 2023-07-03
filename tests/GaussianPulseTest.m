% Define the pulse parameters
A = 1;
x0 = 0;
y0 = 0;
sigma = 1;
k0 = 2*pi;
theta = 0;

% Create the GaussianPulse2D object
pulse = GaussianPulse2D(A, x0, y0, sigma, k0, theta);

% Generate a grid of points
x = linspace(-5, 5, 100);
y = linspace(-5, 5, 100);
[X, Y] = meshgrid(x, y);

% Generate the pulse
P = pulse.generate(X, Y);

% Plot the pulse amplitude
imagesc(x, y, abs(P));
colorbar;
xlabel('x');
ylabel('y');
title('Gaussian pulse amplitude');
