function [X, fval, exitflag, output, grad] = setupMinimization(X0, MinVar, Ns, d, MaxIT)
    % Setup the minimization algorithm
    
    f = @(X) calculateObjective(X, MinVar, Ns, d);
    
    options = optimoptions(@fminunc, 'OutputFcn', @outfun, 'Algorithm', 'quasi-newton', ...
        'Display', 'iter-detailed', 'GradObj', 'on', 'TolFun', 1e-12, ...
        'MaxIter', MaxIT);
    [X, fval, exitflag, output, grad] = fminunc(f, X0, options);
    
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
