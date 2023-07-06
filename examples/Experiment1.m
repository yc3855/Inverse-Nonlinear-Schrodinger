% This experiment reconstructs only k

% What variables am I reconstructing
% MinVar = {'k', 'gamma', 'sigmaTPA'};
MinVar = {'k'};

% Define the discretization of the domain
Nx = 100;
Ny = 100;
x = linspace(0, 1, Nx);
y = linspace(0, 1, Ny);
[X, Y] = meshgrid(x, y);
T = 1;
t = linspace(0, T, Nt);
dx = x(1) - x(0);
dy = y(1) - y(0);
dt = t(1) - t(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% True Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call the setupInitialGuess function
[k_0, gamma_0, sigmaTPA_0, sigma_0] = setupInitialGuess(MinVar, ...
    background_k, k_t, background_gamma, gamma_t, background_sigmaTPA, ...
    sigmaTPA_t, background_sigma, sigma_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Generate sources (initial conditions) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of sources
Ns = 10;

% Create a cell array to store parameters for each source
pulse_params = cell(1, Ns);

% Range for A
A_range = [0.5, 1.5];

% Range for sigma
sigma_range = [0.05, 0.1];

% Iterate to create parameters for each source
for i = 1:Ns
    A = rand * (A_range(2) - A_range(1)) + A_range(1);
    x0 = 0.3 + rand * 0.4; % Random position closer to the center
    y0 = 0.3 + rand * 0.4; % Random position closer to the center
    sigma = rand * (sigma_range(2) - sigma_range(1)) + sigma_range(1);
    k0 = 2 * pi;
    theta = rand * 2 * pi; % Random direction
    pulse_params{i} = [A, x0, y0, sigma, k0, theta];
end

% Generate the sources
sources = generateSources(Ns, pulse_params, X, Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = generateSyntheticData(Nx, Ny, Nt, dx, dy, dt, sources, ...
    k_t, sigma_t, sigmaTPA_t, gamma_t, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Call Minimization %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




