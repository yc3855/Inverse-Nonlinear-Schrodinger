function [wbreal_ret,wbimag_ret]=NLS_adjoint(k,kb,gamma,sigmaTPA,sigma,u_real,u_imag,d_real,d_imag)

% Solve w_t = (i/2k)(w_xx + w_yy) + 2 (i gamma - 1/2 sigma_TPA) |u|^2 w
%            +(-i gamma - 1/2 sigma_TPA) \overline{u}^2 \overline{w} - 1/2 sigma w
% in 2D 0 < x < 1, 0 < y < 1 by DG
% periodic boundary conditions

qx = 1; % the degree for u in x-axis
qy = 1; % the degree for u in y-axis
Nx = 20; % the number of cells in x-axis
Ny = 20; % the number of cells in y-axis
T = 1; % the simulation time
% central flux

% set up the grid
x = linspace(0,1,Nx+1);
hx = x(2:Nx+1) - x(1:Nx);
y = linspace(0,1,Ny+1);
hy = y(2:Ny+1) - y(1:Ny);

% construct the matrices

% get the nodes in reference domain [-1,1] and their weight
[rx,Wx] = GaussQCofs(qx+1);
[ry,Wy] = GaussQCofs(qy+1);

% get the qx+1*qx+1 matrix of legendre polynomials at nodes rx; the maximum degree is q 
Px = zeros(qx+1,qx+1);
for d = 1:qx+1
    Px(:,d) = JacobiP(rx,0,0,d-1);
end

Py = zeros(qy+1,qy+1);
for d = 1:qy+1
    Py(:,d) = JacobiP(ry,0,0,d-1);
end

% get the q+1 * q+1 matrix of the first derivative of legendre polynomials at nodes
% r; the maximum degree is q-1
Sx = zeros(qx+1,qx+1);
for d = 1:qx+1
    Sx(:,d) = GradJacobiP(rx,0,0,d-1);
end

Sy = zeros(qy+1,qy+1);
for d = 1:qy+1
    Sy(:,d) = GradJacobiP(ry,0,0,d-1);
end

% express the boundary [-1;1]
bp = [-1; 1];

Vb_x = zeros(2,qx+1); 
Uxb_x = Vb_x;
for d = 1:qx+1
    Vb_x(:,d) = JacobiP(bp,0,0,d-1);
    Uxb_x(:,d) = GradJacobiP(bp,0,0,d-1); 
end

Vb_y = zeros(2,qy+1); 
Uxb_y = Vb_y;
for d = 1:qy+1
    Vb_y(:,d) = JacobiP(bp,0,0,d-1);
    Uxb_y(:,d) = GradJacobiP(bp,0,0,d-1); 
end

% get the ((qx+1)*(qy+1)) * ((qx+1)*(qy+1)) matrix of legendre polynomials at nodes r; the maximum degree is q 
P = kron(Px,Py);
W = kron(diag(Wx),diag(Wy));

% fluxes matrices

% Initialize w
w_real = zeros((qx+1)*(qy+1),(Nx*Ny));
w_imag = zeros((qx+1)*(qy+1),(Nx*Ny));

% For initial data compute the L2 projections
xloc = zeros(qx+1,Nx);
yloc = zeros(qy+1,Ny);
for i=1:Nx
    xloc(:,i) = (x(i)+x(i+1)+hx(i)*rx)/2;
end

for j=1:Ny
    yloc(:,j) = (y(j)+y(j+1)+hy(j)*ry)/2;
end

% k is an array of shape [(qx+1)*(qy+1), Nx*Ny]
% kb is an array of shape [(qx+1)+(qy+1), Nx*Ny]
% gamma is an array of shape [(qx+1)*(qy+1), Nx*Ny]
% sigmaTPA is an array of shape [(qx+1)*(qy+1), Nx*Ny]
% sigma is an array of shape [(qx+1)*(qy+1), Nx*Ny]

for i=1:Nx
    for j=1:Ny
        
        wloc_real = u_real - d_real;
        wloc_imag =-(u_imag - d_imag);
        
        for d = 1:((qx+1)*(qy+1))
            w_real(d,(i-1)*Ny+j) = (wloc_real*(W))*P(:,d);
            w_imag(d,(i-1)*Ny+j) = (wloc_imag*(W))*P(:,d);
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time stepping - SSPRK(5,4)
CFL = 0.75;
dt = CFL/(2*pi)*(min(hx));
dt = 0.1*dt;
nsteps = ceil(T/dt);
disp(nsteps)
dt = T/nsteps;

wreal_ret = zeros(((qx+1)*(qy+1)),(Nx*Ny),nsteps+1);
wreal_ret(:,:,1) = P*w_real;
wimag_ret = zeros(((qx+1)*(qy+1)),(Nx*Ny),nsteps+1);
wimag_ret(:,:,1) = P*w_imag;

wbreal_ret = zeros((Nx*Ny),nsteps+1);
wbimag_ret = zeros((Nx*Ny),nsteps+1);
Vb = squeeze(kron(Vb_x(1,:), Vb_y(1,:)));
wbreal_ret(:,1) = Vb*w_real;
wbimag_ret(:,1) = Vb*w_imag;

c11 = (0.391752226571890);

c20 = (0.444370493651235);
c21 = (0.555629506348765);
c22 = (0.368410593050371);

c30 = (0.620101851488403);
c32 = (0.379898148511597);
c33 = (0.251891774271694);

c40 = (0.178079954393132);
c43 = (0.821920045606868);
c44 = (0.544974750228521);

c52 = (0.517231671970585);
c53 = (0.096059710526147);
c53_1 = (0.063692468666290);
c54 = (0.386708617503269);
c55 = (0.226007483236906);

for it = 1:nsteps
    
    [rhsw_real, rhsw_imag] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, w_real, ...
        w_imag, k, kb, gamma, sigmaTPA, sigma, P, W, Sx, Sy, Wx, Wy, u_real(:,:,nsteps-it), u_imag(:,:,nsteps-it));
    w_real1 = w_real + c11*dt*rhsw_real;    
    w_imag1 = w_imag + c11*dt*rhsw_imag;
    
    [rhsw_real, rhsw_imag] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, w_real1, ...
        w_imag1, k, kb, gamma, sigmaTPA, sigma, P, W, Sx, Sy, Wx, Wy, u_real(:,:,nsteps-it), u_imag(:,:,nsteps-it));
    w_real2 = c20*w_real + c21*w_real1 + c22*dt*rhsw_real;  
    w_imag2 = c20*w_imag + c21*w_imag1 + c22*dt*rhsw_imag;
    
    
    [rhsw_real, rhsw_imag] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, w_real2, ...
        w_imag2, k, kb, gamma, sigmaTPA, sigma, P, W, Sx, Sy, Wx, Wy, u_real(:,:,nsteps-it), u_imag(:,:,nsteps-it));
    w_real3 = c30*w_real + c32*w_real2 + c33*dt*rhsw_real;   
    w_imag3 = c30*w_imag + c32*w_imag2 + c33*dt*rhsw_imag;
    
    
    [rhsw_real3, rhsw_imag3] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, w_real3, ...
        w_imag3, k, kb, gamma, sigmaTPA, sigma, P, W, Sx, Sy, Wx, Wy, u_real(:,:,nsteps-it), u_imag(:,:,nsteps-it));
    w_real4 = c40*w_real + c43*w_real3 + c44*dt*rhsw_real3;   
    w_imag4 = c40*w_imag + c43*w_imag3 + c44*dt*rhsw_imag3;
    
    
    [rhsw_real, rhsw_imag]   = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, w_real4, ...
        w_imag4, k, kb, gamma, sigmaTPA, sigma, P, W, Sx, Sy, Wx, Wy, u_real(:,:,nsteps-it), u_imag(:,:,nsteps-it));
    w_real  = c52*w_real2 + c53*w_real3 + c53_1*dt*rhsw_real3 + c54*w_real4 + c55*dt*rhsw_real; 
    w_imag  = c52*w_imag2 + c53*w_imag3 + c53_1*dt*rhsw_imag3 + c54*w_imag4 + c55*dt*rhsw_imag;
    
    wreal_ret(:,:,it+1) = P*w_real;
    wimag_ret(:,:,it+1) = P*w_imag;
    
    wbreal_ret(:,it+1) = Vb*w_real;
    wbimag_ret(:,it+1) = Vb*w_imag;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot the projection

plot_wreal = (P*w_real);
plot_wreal_reshape = zeros((qx+1)*(qy+1),Nx*Ny);

plot_wimag = (P*w_imag);
plot_wimag_reshape = zeros((qx+1)*(qy+1),Nx*Ny);

plot_intensity = (P*w_real).^2 + (P*w_imag).^2;
plot_intensity_reshape = zeros((qx+1)*(qy+1),Nx*Ny);

for j=1:Ny
    for i=1:Nx
        plot_wreal_reshape(:,(j-1)*Ny+i) = plot_wreal(:,(i-1)*Ny+j);
        plot_wimag_reshape(:,(j-1)*Ny+i) = plot_wimag(:,(i-1)*Ny+j);
        plot_intensity_reshape(:,(j-1)*Ny+i) = plot_intensity(:,(i-1)*Ny+j);
    end
end

plot_z = zeros((qx+1)*Nx,(qy+1)*Ny);
plot_zimag = zeros((qx+1)*Nx,(qy+1)*Ny);
plot_zintensity = zeros((qx+1)*Nx,(qy+1)*Ny);

for j=1:Ny
    for i=1:Nx
        for k=1:(qy+1)
            new = plot_wreal_reshape(k,(j-1)*Nx+i);
            newimag = plot_wimag_reshape(k,(j-1)*Nx+i);
            newintensity = plot_intensity_reshape(k,(j-1)*Nx+i);
            
            for l=2:(qx+1)
                new(l) = plot_wreal_reshape(k+(l-1)*(qy+1),(j-1)*Nx+i);
                newimag(l) = plot_wimag_reshape(k+(l-1)*(qy+1),(j-1)*Nx+i);
                newintensity(l) = plot_intensity_reshape(k+(l-1)*(qy+1),(j-1)*Nx+i);
            end
            
            plot_z(((i-1)*(qx+1)+1):((i-1)*(qx+1)+qx+1), (j-1)*(qy+1)+k) = new;
            plot_zimag(((i-1)*(qx+1)+1):((i-1)*(qx+1)+qx+1), (j-1)*(qy+1)+k) = newimag;
            plot_zintensity(((i-1)*(qx+1)+1):((i-1)*(qx+1)+qx+1), (j-1)*(qy+1)+k) = newintensity;
        end
    end
end

figure('Position', [50, 50, 400, 1200]);

subplot(3,1,1)
plot_x = reshape(xloc,[],1);
plot_y = reshape(yloc,[],1);
[X, Y] = meshgrid(plot_x, plot_y);
mesh(X, Y, plot_z);
colorbar
title('the real part of the phase')

subplot(3,1,2);
mesh(X, Y, plot_zimag);
colorbar
title('the imaginary part of the phase')

subplot(3,1,3);
mesh(X, Y, plot_zintensity);
% colorbar
title('intensity')

end


function [wt_real, wt_imag] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, w_real, ...
    w_imag, k, kb, gamma, sigmaTPA, sigma, P, W, Px, Py, Sx, Sy, Wx, Wy, u_real, u_imag)

wt_real = zeros((qx+1)*(qy+1),(Nx*Ny));
wt_imag = zeros((qx+1)*(qy+1),(Nx*Ny));

for i=1:Nx
    for j=1:Ny
        
        % local value coressponds to (i,j)
        wreal_local = w_real(:, (i-1)*Ny + j);
        wimag_local = w_imag(:, (i-1)*Ny + j);
        
        ureal_local = squeeze(u_real(:, (i-1)*Ny + j, 1));
        uimag_local = squeeze(u_imag(:, (i-1)*Ny + j, 1));
        
        k_local = k(:, (i-1)*Ny + j);
        gamma_local = gamma(:, (i-1)*Ny + j);
        sigmaTPA_local = sigmaTPA(:, (i-1)*Ny + j);
        sigma_local = sigma(:, (i-1)*Ny + j);
        
        % assemble ut
        Mw1 = (kron(Sx,Py)')*diag(1./(2.*k_local))*W*(kron(Sx,Py));
        Mw1 = Mw1 + (kron(Px,Sy)')*diag(1./(2.*k_local))*W*(kron(Px,Sy));
        
        fu2 = (P*ureal_local).^2 + (P*uimag_local).^2;
        Mw2 = 2*(hx(1)/2)*(hy(1)/2)*(P')*diag(fu2)*diag(gamma_local)*W*P;
        
        Mw3 = (hx(1)/2)*(hy(1)/2)*(P')*diag(fu2)*diag(sigmaTPA_local)*W*P;
        
        fu3_real = (P*ureal_local).^2 - (P*uimag_local).^2;
        fu3_imag = 2*(P*ureal_local).*(P*uimag_local);
        Mw4_real = (hx(1)/2)*(hy(1)/2)*(P')*diag(fu3_real)*diag(gamma_local)*W*P;
        Mw4_imag = (hx(1)/2)*(hy(1)/2)*(P')*diag(fu3_imag)*diag(gamma_local)*W*P;
        
        Mw5_real = (1/2)*(hx(1)/2)*(hy(1)/2)*(P')*diag(fu3_real)*diag(sigmaTPA_local)*W*P;
        Mw5_imag = (1/2)*(hx(1)/2)*(hy(1)/2)*(P')*diag(fu3_imag)*diag(sigmaTPA_local)*W*P;
        
        Mw6 = (1/2)*(hx(1)/2)*(hy(1)/2)*(P')*diag(sigma_local)*W*P;
        
        wt_real(:,(i-1)*Ny+j) = (2/hx(1))*(2/hy(1))*(-Mw1*wimag_local - Mw2*wimag_local ...
            - Mw3*wreal_local - Mw4_real*wimag_local - Mw4_imag*wreal_local ...
            - Mw5_real*wreal_local + Mw5_imag*wimag_local - Mw6*wreal_local);
        
        wt_imag(:,(i-1)*Ny+j) = (2/hx(1))*(2/hy(1))*( Mw1*wreal_local + Mw2*wreal_local ...
            - Mw3*wimag_local - Mw4_real*wreal_local + Mw4_imag*wimag_local ...
            + Mw5_real*wimag_local + Mw5_imag*wreal_local - Mw6*wimag_local);
        
        % central flux: star = (local+out)/2
        % periodic BC
        % (i-1)*Ny + j
        if (j==1) % out corresponds to j == Ny (i-1)*Ny + Ny south side
           wreal_out = w_real(:, (i-1)*Ny + Ny);
           wimag_out = w_imag(:, (i-1)*Ny + Ny);
           kb_s = kb(1:(qx+1), (i-1)*Ny + Ny);
        else % out corresponds to j == j-1 (i-1)*Ny + j - 1
           wreal_out = w_real(:, (i-1)*Ny + j-1);
           wimag_out = w_imag(:, (i-1)*Ny + j-1);
           kb_s = kb(1:(qx+1), (i-1)*Ny + j-1);
        end
        
        Fw_real =         - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Vb_y(1,:)')*Uxb_y(1,:)) * (-wimag_local);
        Fw_real = Fw_real - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Vb_y(1,:)')*Uxb_y(2,:)) * (-wimag_out);
        Fw_real = Fw_real - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Uxb_y(1,:)')*Vb_y(1,:)) * (-wimag_local);
        Fw_real = Fw_real - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Uxb_y(1,:)')*Vb_y(2,:)) * ( wimag_out);
        
        Fw_imag =         - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Vb_y(1,:)')*Uxb_y(1,:)) * ( wreal_local);
        Fw_imag = Fw_imag - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Vb_y(1,:)')*Uxb_y(2,:)) * ( wreal_out);
        Fw_imag = Fw_imag - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Uxb_y(1,:)')*Vb_y(1,:)) * ( wreal_local);
        Fw_imag = Fw_imag - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Uxb_y(1,:)')*Vb_y(2,:)) * (-wreal_out);
        
        if (j==Ny) % out corresponds to j == 1 north side
           wreal_out = w_real(:, (i-1)*Ny + 1);
           wimag_out = w_imag(:, (i-1)*Ny + 1);
           kb_n = kb(1:(qx+1),(i-1)*Ny + 1);
        else % out corresponds to j + 1
           wreal_out = w_real(:, (i-1)*Ny + j+1);
           wimag_out = w_imag(:, (i-1)*Ny + j+1);
           kb_n = kb(1:(qx+1),(i-1)*Ny + j+1);
        end
        
        Fw_real = Fw_real + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Vb_y(2,:)')*Uxb_y(2,:)) * (-wimag_local);
        Fw_real = Fw_real + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Vb_y(2,:)')*Uxb_y(1,:)) * (-wimag_out);
        Fw_real = Fw_real + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Uxb_y(2,:)')*Vb_y(2,:)) * (-wimag_local);
        Fw_real = Fw_real + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Uxb_y(2,:)')*Vb_y(1,:)) * ( wimag_out);
        
        Fw_imag = Fw_imag + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Vb_y(2,:)')*Uxb_y(2,:)) * ( wreal_local);
        Fw_imag = Fw_imag + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Vb_y(2,:)')*Uxb_y(1,:)) * ( wreal_out);
        Fw_imag = Fw_imag + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Uxb_y(2,:)')*Vb_y(2,:)) * ( wreal_local);
        Fw_imag = Fw_imag + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Uxb_y(2,:)')*Vb_y(1,:)) * (-wreal_out);
        
        if (i==1) % out corresponds to i = Nx west side
           wreal_out = w_real(:, (Nx-1)*Ny + j);
           wimag_out = w_imag(:, (Nx-1)*Ny + j);
           kb_w = kb((qx+2):(qx+1+qy+1),(Nx-1)*Ny + j);
        else % out corresponds to i - 1
           wreal_out = w_real(:, (i-2)*Ny + j);
           wimag_out = w_imag(:, (i-2)*Ny + j);
           kb_w = kb((qx+2):(qx+1+qy+1),(i-2)*Ny + j);
        end
        
        Fw_real = Fw_real - (1/2) * kron((Vb_x(1,:)')*Uxb_x(1,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * (-wimag_local);
        Fw_real = Fw_real - (1/2) * kron((Vb_x(1,:)')*Uxb_x(2,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * (-wimag_out);
        Fw_real = Fw_real - (1/2) * kron((Uxb_x(1,:)')*Vb_x(1,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * (-wimag_local);
        Fw_real = Fw_real - (1/2) * kron((Uxb_x(1,:)')*Vb_x(2,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * ( wimag_out);
        
        Fw_imag = Fw_imag - (1/2) * kron((Vb_x(1,:)')*Uxb_x(1,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * ( wreal_local);
        Fw_imag = Fw_imag - (1/2) * kron((Vb_x(1,:)')*Uxb_x(2,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * ( wreal_out);
        Fw_imag = Fw_imag - (1/2) * kron((Uxb_x(1,:)')*Vb_x(1,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * ( wreal_local);
        Fw_imag = Fw_imag - (1/2) * kron((Uxb_x(1,:)')*Vb_x(2,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * (-wreal_out);
        
        if (i==Nx) % out corresponds to i = 1 east side
           wreal_out = w_real(:, j);
           wimag_out = w_imag(:, j);
           kb_e = kb((qx+2):(qx+1+qy+1),j);
        else % out corresponds to i + 1
           wreal_out = w_real(:, i*Ny + j);
           wimag_out = w_imag(:, i*Ny + j);
           kb_e = kb((qx+2):(qx+1+qy+1),i*Ny + j);
        end
        
        Fw_real = Fw_real + (1/2) * kron((Vb_x(2,:)')*Uxb_x(2,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * (-wimag_local);
        Fw_real = Fw_real + (1/2) * kron((Vb_x(2,:)')*Uxb_x(1,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * (-wimag_out);
        Fw_real = Fw_real + (1/2) * kron((Uxb_x(2,:)')*Vb_x(2,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * (-wimag_local);
        Fw_real = Fw_real + (1/2) * kron((Uxb_x(2,:)')*Vb_x(1,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * ( wimag_out);
        
        Fw_imag = Fw_imag + (1/2) * kron((Vb_x(2,:)')*Uxb_x(2,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * ( wreal_local);
        Fw_imag = Fw_imag + (1/2) * kron((Vb_x(2,:)')*Uxb_x(1,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * ( wreal_out);
        Fw_imag = Fw_imag + (1/2) * kron((Uxb_x(2,:)')*Vb_x(2,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * ( wreal_local);
        Fw_imag = Fw_imag + (1/2) * kron((Uxb_x(2,:)')*Vb_x(1,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * (-wreal_out);
        
        Fw_real = Fw_real/2;
        Fw_imag = Fw_imag/2;
        
        wt_real(:,(i-1)*Ny+j) = wt_real(:,(i-1)*Ny+j) + (2/hx(i))*(2/hy(j))*Fw_real;
        wt_imag(:,(i-1)*Ny+j) = wt_imag(:,(i-1)*Ny+j) + (2/hx(i))*(2/hy(j))*Fw_imag;
        
    end
end

end

