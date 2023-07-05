function plotTrueCoefficients(MinVar, k_t, gamma_t, sigmaTPA_t, sigma_t)
    % plotTrueCoefficients Plot the true coefficients k, gamma, sigmaTPA, sigma.
    %
    % Parameters:
    %   MinVar : cell array of strings
    %       The coefficients to be plotted.
    %   k_t : column vector
    %       The true values for the coefficient k.
    %   gamma_t : column vector
    %       The true values for the coefficient gamma.
    %   sigmaTPA_t : column vector
    %       The true values for the coefficient sigmaTPA.
    %   sigma_t : column vector
    %       The true values for the coefficient sigma.

    % Create a figure with 2x2 subplots
    figure;
    sgtitle('True Coefficients');
    
    % Plot k
    subplot(2, 2, 1);
    plot(k_t);
    title('k');
    
    % Plot gamma
    subplot(2, 2, 2);
    plot(gamma_t);
    title('gamma');
    
    % Plot sigmaTPA
    subplot(2, 2, 3);
    plot(sigmaTPA_t);
    title('sigmaTPA');
    
    % Plot sigma
    subplot(2, 2, 4);
    plot(sigma_t);
    title('sigma');
end
