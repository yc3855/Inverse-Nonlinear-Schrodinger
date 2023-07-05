function [X, fval, exitflag, output, grad] = setupMinimization(X0, MinVar, Ns, D, MaxIT)
    % Setup the minimization algorithm
    
    f = @(X) calculateObjective(X, MinVar, Ns, D);
    
    options = optimoptions(@fminunc, 'OutputFcn', @outfun, 'Algorithm', 'quasi-newton', ...
        'Display', 'iter-detailed', 'GradObj', 'on', 'TolFun', 1e-12, ...
        'MaxIter', MaxIT);
    [X, fval, exitflag, output, grad] = fminunc(f, X0, options);
    
end
