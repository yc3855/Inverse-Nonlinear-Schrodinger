function plotFinalResults(MinVar, X)
    % Reconstructed coefficients
    k_r = X(1:M);
    gamma_r = X(M + 1 : 2 * M);
    sigmaTPA_r = X(2 * M + 1 : 3 * M);
    sigma_r = X(3 * M + 1 : 4 * MinVar);
    
    % Plot reconstruction results
    if ismember("k", MinVar)
        % plot k_r
    end
    if ismember("gamma", MinVar)
        % plot gamma_r
    end
    if ismember("sigmaTPA", MinVar)
        % plot sigmaTPA_r
    end
    if ismember("sigma", MinVar)
        % plot sigma_r
    end
end