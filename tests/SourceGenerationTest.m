% Define the domain
x = linspace(0, 1, 100);
y = linspace(0, 1, 100);
[X, Y] = meshgrid(x, y);

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

% Visualize the real part of the sources
for i = 1:Ns
    subplot(ceil(sqrt(Ns)), ceil(sqrt(Ns)), i);
    imagesc(x, y, real(sources{i}));
    title(['Source ', num2str(i)]);
    colorbar;
    axis equal;
    xlim([0, 1]);
    ylim([0, 1]);
end
