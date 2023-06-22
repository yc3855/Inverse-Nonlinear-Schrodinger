% Coefficients to reconstruct; subset of ["k","gamma","sigmaTPA", "sigma"]
MinVar=["gamma"]; 

noiselevel=0.00; % set noise level
beta_k=0e-7; beta_gamma=0e-7; beta_sigmaTPA=0e-7; beta_sigma=0; regularization parameters

MaxIT=1000; % max number of iterations in optimization algorithm

%iter_plot=50; % plots intermediate reconstructions after every iter_plot iterations

Ns=36; % number of sources (initial conditions)

%Set true coefficients
%True k
k_t = Profile(); 
k_t.background = 1;
k_t.rectangles = [Rectangle([0.4 1.3; 1 1.5],0.5)];

%True gamma
gamma_t = Profile(); 
gamma_t.background = 0.1;
gamma_t.rectangles = [Rectangle([0.4 0.8; 0.2 1.1],0.03), ...
    Rectangle([0.6 1.5; 1.5 1.8],0.06), Rectangle([1.3 1.7; 0.4 0.8],0.1)];

%True sigmaTPA
sigmaTPA_t = Profile();
sigmaTPA_t.background = 0.2;
sigmaTPA_t.rectangles = [Rectangle([0.5 0.8; 0.5 0.8],0.2), ...
    Rectangle([1.4 1.7; 1.4 1.7],0.4), Rectangle([1.0 1.8; 0.3 0.6],0.6)];

%True sigma
sigma_t = Profile();
sigma_t.background = 1;
sigma_t.circles = [Circle([1.5 1.0],0.3,1)];
sigma_t.rectangles = [Rectangle([0.5 0.8; 0.5 0.8],2)];


% Load geometrical information on the domain
load geo-2b2

NLS(MinVar, k_t, gamma_t, sigmaTPA_t, sigma_t, Ns, noiselevel, MaxIT, ...
    beta_k, beta_gamma, beta_sigmaTPA, beta_sigma)