function X0 = setupInitialGuess(MinVar, k_t, gamma_t, sigmaTPA_t, sigma_t)
    % Set up initial guesses
    % k_0 = background of k_t
    % gamma_0 = background of gamma_t
    % sigmaTPA_0 = background of sigmaTPA_t
    % sigma_0 = background of sigma_t

    % Set up initial guesses only for coefficients we want to reconstruct
    if ~ismember("k",MinVar)
        k_0 = k_t;
    end

    if ~ismember("gamma",MinVar)
        gamma_0 = gamma_t;
    end

    if ~ismember("sigmaTPA",MinVar)
        sigmaTPA_0 = sigmaTPA_t;
    end

    if ~ismember("sigma",MinVar)
        sigma_0 = sigma_t;
    end

    X0=[k_0' gamma_0' sigmaTPA_0' sigma_0']';
end
