function d = generateSyntheticData(M, Ns, noiselevel, T, F, k_t, gamma_t, sigmaTPA_t, sigma_t)
    % generateSyntheticData Generate noisy data from the solution of the nonlinear
    % Schrödinger equation.
    %
    % Parameters:
    %   M : int
    %       The number of measurements.
    %   Ns : int
    %       The number of sources.
    %   noiselevel : float
    %       The level of noise to add to the data.
    %   T : float
    %       The final time at which we want the solution.
    %   F : 2D array (M x Ns)
    %       The Ns sources/initial conditions for the NLS equation
    %   k_t : float
    %       Coefficient for the NLS equation.
    %   gamma_t : float
    %       Coefficient for the NLS equation.
    %   sigmaTPA_t : float
    %       Coefficient for the NLS equation.
    %   sigma_t : float
    %       Coefficient for the NLS equation.
    %
    % Returns:
    %   d : 2D array (M x Ns)
    %       The generated data, with added noise.

    % Preallocate the data array
    d = zeros(M, Ns);

    % Loop over the number of sources
    for s = 1:Ns
        f_s = F(:,s);
        % Run NLS equation with initial condition f_s and coefficients
        d_s = NLS_forward(T, f_s, k_t, gamma_t, sigmaTPA_t, sigma_t);

        % Add noise to data
        d(:, s) = d_s * (1 + noiselevel * 2 * (rand(M, 1) - 0.5));
    end
end
