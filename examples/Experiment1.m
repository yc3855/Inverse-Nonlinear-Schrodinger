% Define the discretization of the domain
x = linspace(0, 1, 100);
y = linspace(0, 1, 100);

% Reflective Index (k)
background_k = 1.5; % Typical for glass
rectangles_k = [Rectangle([0.1, 0.1], [0.2, 0.2], 1.33)]; % Typical for water
circles_k = [Circle([0.7, 0.7], 0.1, 1.0)]; % Typical for air

% Kerr Coefficient (gamma)
background_gamma = 1e-20; % Typical for glass in m²/W
rectangles_gamma = [Rectangle([0.3, 0.4], [0.4, 0.5], 5e-20)];
circles_gamma = [Circle([0.6, 0.6], 0.08, 3e-19)];

% Rescaled Two-Photon Absorption Coefficient (sigmaTPA)
background_sigmaTPA = 1e-5; % in cm/W
rectangles_sigmaTPA = [Rectangle([0.5, 0.6], [0.6, 0.7], 2e-5)];
circles_sigmaTPA = [Circle([0.4, 0.4], 0.07, 3e-5)];

% Linear Physical Absorption Coefficient (sigma)
background_sigma = 100; % in m⁻¹
rectangles_sigma = [Rectangle([0.7, 0.7], [0.8, 0.8], 200)];
circles_sigma = [Circle([0.5, 0.5], 0.1, 150)];

% Generate k parameter
k_t = generate_k(x, y, background_k, rectangles_k, circles_k);

% Generate sigma parameter
sigma_t = generate_sigma(x, y, background_sigma, rectangles_sigma, circles_sigma);

% Generate sigmaTPA parameter
sigmaTPA_t = generate_sigmaTPA(x, y, background_sigmaTPA, rectangles_sigmaTPA, circles_sigmaTPA);

% Generate sigmaTPA parameter
gamma_t = generate_gamma(x, y, background_gamma, rectangles_gamma, circles_gamma);

figure;

% Plot for k_t parameter
subplot(2, 2, 1); % This creates a 2-row, 2-column grid of plots and selects the first one for plotting.
imagesc(x, y, k_t);
colorbar;
title('k parameter');
axis equal;

% Plot for gamma_t parameter
subplot(2, 2, 2); % This selects the second plot for plotting.
imagesc(x, y, gamma_t);
colorbar;
title('gamma parameter');
axis equal;

% Plot for sigmaTPA_t parameter
subplot(2, 2, 3); % This selects the third plot for plotting.
imagesc(x, y, sigmaTPA_t);
colorbar;
title('sigmaTPA parameter');
axis equal;

% Plot for sigma_t parameter
subplot(2, 2, 4); % This selects the fourth plot for plotting.
imagesc(x, y, sigma_t);
colorbar;
title('sigma parameter');
axis equal;
