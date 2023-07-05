function plotInitialGuess(MinVar, X0)
    % plotInitialGuess Plot the initial guesses for the coefficients k, gamma, sigmaTPA, sigma.
    %
    % Parameters:
    %   MinVar : cell array of strings
    %       The coefficients to be plotted.
    %   X0 : column vector
    %       The initial guesses for the coefficients.

    % Reshape the initial guess vector into 2D
    Nx = sqrt(length(X0));
    Ny = Nx;
    X0 = reshape(X0, Nx, Ny);

    % Create a figure with 2x2 subplots
    figure;
    sgtitle('The Initial Guess');
    for i = 1:length(MinVar)
        subplot(2, 2, i);
        parameter = MinVar{i};
        plot(X0);
        title(parameter, 'Interpreter', 'latex');
    end
end
