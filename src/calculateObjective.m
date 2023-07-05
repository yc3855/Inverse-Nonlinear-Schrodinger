function [Phi, grad] = calculateObjective(X, MinVar, Ns, d)
    % NLSObj Evaluate the objective function and its gradients.
    %
    % Parameters:
    %   X : column vector
    %       The optimization variables (k, gamma, sigmaTPA, sigma).
    %   MinVar : cell array of strings
    %       The coefficients to be optimized.
    %   Ns : int
    %       The number of sources.
    %   d : 2D array (M x Ns)
    %       The observed data.

    M = Nx * Ny; % Total number of nodes in the spacial mesh

    % Extract the current values of the coefficients
    k_c = X(1:M);
    gamma_c = X(M+1:2*M);
    sigmaTPA_c = X(2*M+1:3*M);
    sigma_c = X(3*M+1:4*M);

    % Initialize the objective function and gradient
    Phi = 0.0;
    grad = zeros(4*M, 1);

    for s = 1:Ns
        % Run NLS equation with current coefficients and forward solver
        [us_real, us_imag, dt] = NLS_forward(k_c, gamma_c, sigmaTPA_c, sigma_c, F(:,s), T);
        us = us_real + 1i*us_imag;
        ds_real = us_real(:, end);
        ds_imag = us_imag(:, end);
        ds = ds_real + 1i * ds_imag;

        rz = ds - d(:, s); % Residual on mesh point locations

        % Contribution to the objective function from source s
        Phi = Phi + 0.5 * sum(abs(rz).^2) * dx * dy;

        if nargout > 1
            % Solve the adjoint equations and compute the gradients

            % Solution to the adjoint terminal value problem
            [ws_real, ws_imag, dt] = NLS_adjoint(k_c, gamma_c, sigmaTPA_c, sigma_c, ...
                us_real, us_imag, real(d(:,s)), imag(d(:,s)), T);
            ws = ws_real + 1i * ws_imag;
            ws = flip(ws, 2);

            % Compute the gradients for each coefficient
            if ismember('k', MinVar)
                grad(1:M) = grad(1:M) + computeGradient('k', k_c, us, ws, dt);
            end

            if ismember('gamma', MinVar)
                grad(M+1:2*M) = grad(M+1:2*M) + computeGradient('gamma', gamma_c, us, ws, dt);
            end

            if ismember('sigmaTPA', MinVar)
                grad(2*M+1:3*M) = grad(2*M+1:3*M) + computeGradient('sigmaTPA', sigmaTPA_c, us, ws, dt);
            end

            if ismember('sigma', MinVar)
                grad(3*M+1:4*M) = grad(3*M+1:4*M) + computeGradient('sigma', sigma_c, us, ws, dt);
            end
        end
    end
end

function gradient = computeGradient(coeff, coeff_c, us, ws, dt)
    % Compute the gradient for a specific coefficient
    
    gradient = sum(real(computeTerm(coeff, coeff_c, us, ws)),2) * dt * dx * dy;
end

function term = computeTerm(coeff, coeff_c, us, ws)
    % Compute a specific term for the gradient calculation
    
    switch coeff
        case 'k'
            % Compute the Laplacian
            u_2D = reshape(us_real, [Nx, Ny]);
            u_Lap = del2(u_2D, dx);
            u_Lap_1D = reshape(u_Lap, [M, 1]);
            
            term = -1i ./ (2 * coeff_c.^2) .* u_Lap_1D .* ws;
        case 'gamma'
            term = 1i * abs(us).^2 .* us .* ws;
        case 'sigmaTPA'
            term = -1/2 * abs(us).^2 .* us .* ws;
        case 'sigma'
            term = -1/2 * us .* ws;
        otherwise
            error('Invalid coefficient specified.');
    end
end
