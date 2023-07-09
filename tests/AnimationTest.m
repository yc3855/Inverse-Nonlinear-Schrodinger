% Set the dimensions and number of time steps
Nx = 100;   % Number of points along the x-axis
Ny = 100;   % Number of points along the y-axis
Nt = 100;   % Number of time steps

% Generate random u_real and u_imag data
u_real = rand(Nx, Ny, Nt);
u_imag = rand(Nx, Ny, Nt);

% Define the step size for the animation
step = 10;  % Plot every 10th frame

% Plot the wavefield animation
plotWavefieldAnimation(u_real, u_imag, step, 'amplitude');
