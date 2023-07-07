function [k, gamma, sigmaTPA, sigma, fval, exitflag, output, grad] = setupMinimization(k_0, ...
    gamma_0, sigmaTPA_0, sigma_0, MinVar, Ns, D, MaxIT)
    % Setup the minimization algorithm
    
    % Matlab optimizer expects a column vector for which to be minimized
    k_0_flat = reshape(k_0, [], 1);  % Flatten k_0 into a column vector
    gamma_0_flat = reshape(gamma_0, [], 1);  % Flatten gamma_0 into a column vector
    sigmaTPA_0_flat = reshape(sigmaTPA_0, [], 1);  % Flatten sigmaTPA_0 into a column vector
    sigma_0_flat = reshape(sigma_0, [], 1);  % Flatten sigma_0 into a column vector

    X0 = [k_0_flat; gamma_0_flat; sigmaTPA_0_flat; sigma_0_flat];  % Concatenate into X0

    % Nx, Ny, M should be in the workspace
    % Determine the dimensions Nx and Ny based on the size of k_0
    [Nx, Ny] = size(k_0);
    M = Nx * Ny; % Total number of nodes in the spacial mesh

    % Lambda function: f(X) = calculateObjective(X, MinVar, Ns, D)
    f = @(X) calculateObjective(X, MinVar, Ns, D, M, sources);
    
    % Set up optimizer
    options = optimoptions(@fminunc, 'OutputFcn', @outfun, 'Algorithm', 'quasi-newton', ...
        'Display', 'iter-detailed', 'GradObj', 'on', 'TolFun', 1e-12, ...
        'MaxIter', MaxIT);
    
    [X, fval, exitflag, output, grad] = fminunc(f, X0, options);
    
    % Properly reshape back to 2D array
    k = reshape(X(1:Nx*Ny), Nx, Ny);  % Reshape the first part of X to match k_0 shape
    gamma = reshape(X(Nx*Ny+1:2*Nx*Ny), Nx, Ny);  % Reshape the second part of X to match gamma_0 shape
    sigmaTPA = reshape(X(2*Nx*Ny+1:3*Nx*Ny), Nx, Ny);  % Reshape the third part of X to match sigmaTPA_0 shape
    sigma = reshape(X(3*Nx*Ny+1:end), Nx, Ny);  % Reshape the fourth part of X to match sigma_0 shape

end
