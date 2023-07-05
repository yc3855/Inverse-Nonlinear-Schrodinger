function NLS(MinVar, k_t, gamma_t, sigmaTPA_t, sigma_t, Ns, noiselevel, MaxIT, ...
    beta_k, beta_gamma, beta_sigmaTPA, beta_sigma)

tic; tb=toc;


M=Nx*Ny; % total number of nodes in the spacial mesh


% Set up initial guesses
% k_0 = background of k_t
% gamma_0 = background of gamma_t
% sigmaTPA_0 = background of sigmaTPA_t
% sigma_0 = background of sigma_t


% Set up initial guesses only for coefficients we want to reconstruct
if ~ismember("k",MinVar)
    k_0 = k_t;
end

if ~ismember("gamma",MinVar)
    gamma_0 = gamma_t;
end

if ~ismember("sigmaTPA",MinVar)
    sigmaTPA_0 = sigmaTPA_t;
end

if ~ismember("sigma",MinVar)
    sigma_0 = sigma_t;
end

X0=[k_0' gamma_0' sigmaTPA_0' sigma_0']';



%Plot initial guesses
if ismember("k",MinVar)
    % plot k_0
end
if ismember("gamma",MinVar)
    % plot gamma_0
end
if ismember("sigmaTPA",MinVar)
    % plot sigmaTPA_0
end
if ismember("sigma",MinVar)
    % plot sigma_0
end



%Plot true coefficients
if ismember("k",MinVar)
    % plot k_t
end
if ismember("gamma",MinVar)
    % plot gamma_t
end
if ismember("sigmaTPA",MinVar)
    % plot sigmaTPA_t
end
if ismember("sigma",MinVar)
    % plot sigma_t
end


% Generating synthetic data
disp(' ');
disp(' ');
disp('Generating synthetic data .......');
disp(' ');

d = GenerateData(M, Ns, noiselevel, T, f_s, k_t, gamma_t, sigmaTPA_t, sigma_t);

disp('Finished generating synthetic data .......');



% Setup the minimization algorithm
disp(' ');
disp(' ');
disp('Minimizing objective function .......');
disp(' ');



f = @(X) NLSObj(X,MinVar,Ns,d);

options = optimoptions(@fminunc,'OutputFcn',@outfun,'Algorithm','quasi-newton', ...
    'Display','iter-detailed','GradObj','on','TolFun',1e-12,...
    'MaxIter',MaxIT);
[X,fval,exitflag,output,grad] = fminunc(f,X0,options);


disp(' ');
disp(' ');
disp('Finished minimizing objective function .......');

disp(' ');
disp(' ');
disp('Plotting final results .......');
disp(' ');

% Reconsructed coefficients
k_r = X(1:M);
gamma_r = X(M+1:2*M);
sigmaTPA_r = X(2*M+1:3*M);
sigma_r = X(3*M+1:4*MinVar);

% Plot reconstruction results
if ismember("k",MinVar)
    % plot k_r
end
if ismember("gamma",MinVar)
    % plot gamma_r
end
if ismember("sigmaTPA",MinVar)
    % plot sigmaTPA_r
end
if ismember("sigma",MinVar)
    % plot sigma_r
end



disp('Finished plotting final results .......');

save Exp01-Info geo P E T SrcInfo BdaryInfo wnum Ns MaxIT ...
                  OptimMethod noiselevel dx dy -ASCII
save Exp01-Results k_t k_0 k_r gamma_t gamma_0 gamma_r ...
    sigmaTPA_t sigmaTPA_0 sigmaTPA_r sigma_t sigma_0 sigma_r -ASCII

te=toc;
disp(' ');
disp(' ');
disp(['The code ran for: ' num2str(te-tb) ' seconds']);
disp(' ');
disp(' ');

end