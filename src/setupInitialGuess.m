function [k_0, gamma_0, sigmaTPA_0, sigma_0] = setupInitialGuess(MinVar, background_k, k_t, background_gamma, gamma_t, background_sigmaTPA, sigmaTPA_t, background_sigma, sigma_t)
    % Set up initial guesses
    % If a coefficient is in MinVar, its initial guess is set to its background value.
    % If a coefficient is not in MinVar, its initial guess is set to its true profile.

    % Initialize the initial guesses with true profiles
    k_0 = k_t;
    gamma_0 = gamma_t;
    sigmaTPA_0 = sigmaTPA_t;
    sigma_0 = sigma_t;

    % Set initial guesses to background values for coefficients in MinVar
    if ismember("k", MinVar)
        k_0 = background_k;
    end

    if ismember("gamma", MinVar)
        gamma_0 = background_gamma;
    end

    if ismember("sigmaTPA", MinVar)
        sigmaTPA_0 = background_sigmaTPA;
    end

    if ismember("sigma", MinVar)
        sigma_0 = background_sigma;
    end
end
