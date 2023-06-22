%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluate the objective function and its gradients with 
% respect to the optimization variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f g] = NLSObj(X,MinVar,Ns,d)

M=Nx*Ny; % total number of nodes in the spacial mesh
% ne = size(SrcInfo,2); % number of edges/nodes on the domain boundary

k_c = X(1:M);% current value of k
gamma_c = X(M+1:2*M); % current value of gamma
sigmaTPA_c = X(2*M+1:3*M); %current value of sigmaTPA
sigma_c = X(3*M+1:4*M); %current value of sigma

f = 0.0;
g = zeros(4*M,1);

for s = 1:Ns
 
    % Run NLS equation with initial condition f_s with coefficients 
    % k_c, gamma_c, sigmaTPA_c, sigma_c using forward solver
    %
    % Let u_s denote the solution for all t in [0,T]
    %
    % Let d_s be the solution at time T

    rz = zeros(M,1); % residual on mesh point locations
    rz = d_s - d(:,s);
    
    % the contribution to the objective function from source s
    f = f + 0.5*sum(abs(rz).^2)*dx*dy;
    
    % the contribution to the gradient from source s
    if nargout > 1         

        % solve the adjoint equations
        %
        % Let w_s be the solution to the adjoint terminal value problem
        % (see overleaf document), using the forward solver for the linear 
        % Schrodinger equation. The terminal condition at time T is conj(rz)

        % the gradient w.r.t k            
        if ismember("k",MinVar)
            %g(1:M) = g(1:M) + [integral_0^T real(-1i./(2*k.^2).*Delta u_s.* w_s) dt]*dx*dy;
            % use naive quadrature rule for the time integral
        end
        % the gradient w.r.t gamma            
        if ismember("gamma",MinVar)
            %g(M+1:2*M) = g(M+1:2*M) + [integral_0^T real(1i*abs(u_s).^2.*u_s.*w_s) dt]*dx*dy;
            % use naive quadrature rule for the time integral
        end
        % the gradient w.r.t sigmaTPA            
        if ismember("sigmaTPA",MinVar)
            %g(2*M+1:3*M) = g(2*M+1:3*M) + [integral_0^T real(-1/2*abs(u_s).^2.*u_s.*w_s) dt]*dx*dy;
            % use naive quadrature rule for the time integral
        end
        % the gradient w.r.t sigma            
        if ismember("sigma",MinVar)
            %g(3*M+1:4*M) = g(3*M+1:4*M) + [integral_0^T real(-1/2*u_s.*w_s) dt]*dx*dy;
            % use naive quadrature rule for the time integral
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