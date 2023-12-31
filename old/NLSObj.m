%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluate the objective function and its gradients with 
% respect to the optimization variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi grad] = NLSObj(X,MinVar,Ns,d)

M=Nx*Ny; % total number of nodes in the spacial mesh
% ne = size(SrcInfo,2); % number of edges/nodes on the domain boundary

k_c = X(1:M);% current value of k
gamma_c = X(M+1:2*M); % current value of gamma
sigmaTPA_c = X(2*M+1:3*M); %current value of sigmaTPA
sigma_c = X(3*M+1:4*M); %current value of sigma

Phi = 0.0;
grad = zeros(4*M,1);

for s = 1:Ns
 
    % Run NLS equation with initial condition f_s with coefficients 
    % k_c, gamma_c, sigmaTPA_c, sigma_c using forward solver
    %
    % Let u_s denote the solution for all t in [0,T]
    %
    % Let d_s be the solution at time T
    [us_real,us_imag,dt]=NLS_forward(k_c,gamma_c,sigmaTPA_c,sigma_c,F(:,s),T);
    us = us_real + 1i*us_imag;
    ds_real = us_real(:,end);
    ds_imag = us_imag(:,end);
    ds = ds_real + 1i*ds_imag;

    rz = zeros(M,1); % residual on mesh point locations
    rz = ds - d(:,s);
    
    % the contribution to the objective function from source s
    Phi = Phi + 0.5*sum(abs(rz).^2)*dx*dy;
    
    % the contribution to the gradient from source s
    if nargout > 1         

        % solve the adjoint equations
        %
        % Let w_s be the solution to the adjoint terminal value problem
        % (see overleaf document), using the forward solver for the linear 
        % Schrodinger equation. For the coefficients, use k_c, gamma_c, ...
        % sigmaTPA_c, sigma_c, and u_s. The terminal condition ...
        % at time T is conj(rz)
        [ws_real,ws_imag,dt]=NLS_adjoint(k_c,gamma_c,sigmaTPA_c,sigma_c,...
            us_real,us_imag,real(d(:,s)),imag(d(:,s)),T);
        ws = ws_real + 1i*ws_imag;
        ws = flip(ws,2);
        

        % use naive rectangle quadrature rule for the time integrals
        % the gradient w.r.t k            
        if ismember("k",MinVar)
            % Assume u is your 1D vector of size Nx*Ny
            u_2D = reshape(u, [Nx, Ny]); % reshape it into a 2D matrix of size Nx by Ny
            
            % Compute the Laplacian
            u_Lap = del2(u_2D, dx);
            
            % Reshape the 2D matrix back into a 1D vector
            u_Lap_1D = reshape(u_Lap, [Nx*Ny, 1]);

            grad(1:M) = grad(1:M) + sum(real(-1i./(2*k.^2).*u_Lap_1D.* w_s), 2)*dt*dx*dy;
        end
        % the gradient w.r.t gamma            
        if ismember("gamma",MinVar)
            grad(M+1:2*M) = grad(M+1:2*M) + sum(real(1i*abs(us).^2.*us.*ws),2)*dt*dx*dy;
        end
        % the gradient w.r.t sigmaTPA            
        if ismember("sigmaTPA",MinVar)
            grad(2*M+1:3*M) = grad(2*M+1:3*M) + sum(real(-1/2*abs(u_s).^2.*u_s.*w_s),2)*dt*dx*dy;
        end
        % the gradient w.r.t sigma            
        if ismember("sigma",MinVar)
            grad(3*M+1:4*M) = grad(3*M+1:4*M) + sum(real(-1/2*u_s.*w_s),2)*dt*dx*dy;
        end
        
    end
    
end




% Add regularization terms to both the objective function and its gradients

% Ignore below, will implement later

% if ismember("Ref", MinVar)
%     [Rx,Ry] = pdegrad(P,T,refc);
%     Rx1=pdeprtni(P,T,Rx); Ry1=pdeprtni(P,T,Ry);
%     f=f+0.5*betan*sum(Rx1.^2+Ry1.^2)*dx*dy;
%     if nargout > 1
%         [Rxx, Rxy]=pdegrad(P,T,Rx1); [Ryx, Ryy]=pdegrad(P,T,Ry1);
%         Rx2=pdeprtni(P,T,Rxx); Ry2=pdeprtni(P,T,Ryy);
%         Deltan=Rx2+Ry2;
%         g(1:M)=g(1:M)-betan*Deltan*dx*dy;
%         for j=1:ne
%             nd=BdaryInfo(1,j);
%             g(nd)=g(nd)+betan*(BdaryInfo(3,j)*Rx1(nd)+BdaryInfo(4,j)*Ry1(nd))*BdaryInfo(5,j);
%         end
%     end
% end
% if ismember("Sigma", MinVar)
%     [Sx,Sy] = pdegrad(P,T,sigmac);
%     Sx1=pdeprtni(P,T,Sx); Sy1=pdeprtni(P,T,Sy);
%     f=f+0.5*betaS*sum(Sx1.^2+Sy1.^2)*dx*dy;
%     if nargout > 1
%         [Sxx, Sxy]=pdegrad(P,T,Sx1); [Syx, Syy]=pdegrad(P,T,Sy1);
%         Sx2=pdeprtni(P,T,Sxx); Sy2=pdeprtni(P,T,Syy);
%         DeltaSigma=Sx2+Sy2;
%         g(M+1:2*M)=g(M+1:2*M)-betaS*DeltaSigma*dx*dy;
%         for j=1:ne
%             nd=BdaryInfo(1,j);
%             g(M+nd)=g(M+nd)+betaS*(BdaryInfo(3,j)*Sx1(nd)+BdaryInfo(4,j)*Sy1(nd))*BdaryInfo(5,j);
%         end
%     end
% end
% if ismember("gamma", MinVar)
%     [gx,gy] = pdegrad(P,T,gammac);
%     gx1=pdeprtni(P,T,gx); gy1=pdeprtni(P,T,gy);
%     f=f+0.5*betag*sum(gx1.^2+gy1.^2)*dx*dy;
%     if nargout > 1
%         [gxx, gxy]=pdegrad(P,T,gx1); [gyx, gyy]=pdegrad(P,T,gy1);
%         gx2=pdeprtni(P,T,gxx); gy2=pdeprtni(P,T,gyy);
%         Deltagamma=gx2+gy2;
%         g(2*M+1:3*M)=g(2*M+1:3*M)-betag*Deltagamma*dx*dy;
%         for j=1:ne
%             nd=BdaryInfo(1,j);
%             g(2*M+nd)=g(2*M+nd)+betag*(BdaryInfo(3,j)*gx1(nd)+BdaryInfo(4,j)*gy1(nd))*BdaryInfo(5,j);
%         end
%     end
% end