function sources = generateSources(Ns, pulse_params, X, Y)
    % Generate a cell array of sources
    % Ns: Number of sources
    % pulse_params: A cell array containing parameters [A, x0, y0, sigma, k0, theta] for each source
    % X, Y: Meshgrid for the domain

    sources = cell(1, Ns);

    for i = 1:Ns
        % Get parameters for the Gaussian pulse
        params = pulse_params{i};
        A = params(1);
        x0 = params(2);
        y0 = params(3);
        sigma = params(4);
        k0 = params(5);
        theta = params(6);

        % Create the GaussianPulse2D object
        pulse = GaussianPulse2D(A, x0, y0, sigma, k0, theta);

        % Generate the pulse
        sources{i} = pulse.generate(X, Y);
    end
end
