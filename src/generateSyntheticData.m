function d = generateSyntheticData(Nx, Ny, Nt, dx, dy, dt, sources, ...
    k_t, sigma_t, sigmaTPA_t, gamma_t, noise)
    % generateSyntheticData Generate noisy data from the solution of the nonlinear
    % Schrödinger equation.
    %
    % Parameters:
    %   Nx, Ny : int
    %       The number of spatial discretization points along the x and y axes.
    %   Nt : int
    %       The number of temporal discretization points.
    %   dx, dy, dt : float
    %       The step sizes for spatial and temporal discretization.
    %   sources : cell array
    %       The cell array containing the source functions.
    %   k_t, gamma_t, sigmaTPA_t, sigma_t : float
    %       The true parameters for the NLS equation.
    %   noise : float
    %       The level of noise to add to the data.
    %
    % Returns:
    %   d : 4D array (Nx x Ny x Nt x Ns)
    %       The generated data, with added noise.

    % Get the number of sources
    Ns = length(sources);

    % Preallocate the data array
    d = zeros(Nx, Ny, Nt, Ns);

    % Loop over the number of sources
    for s = 1:Ns
        source_s = sources{s};
        % Run NLS equation with initial condition source_s and coefficients
        d_s = forward_NLS(Nx, Ny, Nt, dx, dy, dt, source_s, k_t, gamma_t, sigmaTPA_t, sigma_t);

        % Add noise to data
        d(:, :, :, s) = d_s .* (1 + noise * 2 * (rand(Nx, Ny, Nt) - 0.5));
    end
end