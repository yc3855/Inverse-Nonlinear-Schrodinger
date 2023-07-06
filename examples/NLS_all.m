function NLS_all(MinVar, k_t, gamma_t, sigmaTPA_t, sigma_t, Ns, noiselevel, MaxIT, ...
    beta_k, beta_gamma, beta_sigmaTPA, beta_sigma)

tic; tb=toc;
M=Nx*Ny; % total number of nodes in the spacial mesh

% Set up initial guesses
X0 = setupInitialGuess(MinVar, k_t, gamma_t, sigmaTPA_t, sigma_t);

% Generate sources (initial conditions)
F = generateSources(Ns);

% Generating synthetic data
D = generateSyntheticData(M, Ns, noiselevel, T, F, k_t, gamma_t, sigmaTPA_t, sigma_t);

% Setup the minimization algorithm
[X, fval, exitflag, output, grad] = setupMinimization(X0, MinVar, Ns, D, MaxIT);

% Save results
saveResults(MinVar, X, k_t, gamma_t, sigmaTPA_t, sigma_t, X0, geo, P, E, T, SrcInfo, BdaryInfo, wnum, Ns, MaxIT, ...
                  OptimMethod, noiselevel, dx, dy);

te=toc;
disp(['The code ran for: ' num2str(te-tb) ' seconds']);

end
