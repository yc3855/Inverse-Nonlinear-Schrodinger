function [ureal_ret,uimag_ret]=NLSINV_forward(k,kb,gamma,sigmaTPA,sigma)

% Solve u_t = (i/2k(x))(u_xx + u_yy) + (i gamma(x) - 1/2 sigma_TPA(x)) |u|^2 u
%            -1/2 sigma(x) u
% in 2D 0 < x < 1, 0 < y < 1 by DG
% periodic boundary conditions

qx = 1; % qx = degree for u in x-axis
qy = 1; % qy = degree for u in y-axis
Nx = 20; % Nx = number of cells in x-axis
Ny = 20; % Ny = number of cells in y-axis
T = 1; % T is the simulation time

% set up the grid
x = linspace(0,1,Nx+1);
hx = x(2:Nx+1) - x(1:Nx);
y = linspace(0,1,Ny+1);
hy = y(2:Ny+1) - y(1:Ny);

% construct the matrices

% get the nodes in reference domain [-1,1] and their weight
[rx,Wx] = GaussQCofs(qx+1);
[ry,Wy] = GaussQCofs(qy+1);

% get the (qx+1)*(qx+1) matrix of legendre polynomials at nodes rx; the maximum degree is qx 
Px = zeros(qx+1,qx+1);
for d = 1:qx+1
    Px(:,d) = JacobiP(rx,0,0,d-1);
end

% get the (qy+1)*(qy+1) matrix of legendre polynomials at nodes ry; the maximum degree is qy
Py = zeros(qy+1,qy+1);
for d = 1:qy+1
    Py(:,d) = JacobiP(ry,0,0,d-1);
end

% get the (qx+1)*(qx+1) matrix of the first derivative of legendre polynomials at nodes
% rx; the maximum degree is qx-1
Sx = zeros(qx+1,qx+1);
for d = 1:qx+1
    Sx(:,d) = GradJacobiP(rx,0,0,d-1);
end

% get the (qy+1)*(qy+1) matrix of the first derivative of legendre polynomials at nodes
% ry; the maximum degree is qy-1
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

% Initialize u and v
u_real = zeros((qx+1)*(qy+1),(Nx*Ny));
u_imag = zeros((qx+1)*(qy+1),(Nx*Ny));

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
        
        [uloc] = solnew(xloc(:,i),yloc(:,j),0);

        for d = 1:((qx+1)*(qy+1))
            u_real(d,(i-1)*Ny+j) = (uloc(1,:)*(W))*P(:,d);
            u_imag(d,(i-1)*Ny+j) = (uloc(2,:)*(W))*P(:,d);
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

ureal_ret = zeros(((qx+1)*(qy+1)),(Nx*Ny),nsteps+1);
ureal_ret(:,:,0) = P*u_real;
uimag_ret = zeros(((qx+1)*(qy+1)),(Nx*Ny),nsteps+1);
uimag_ret(:,:,0) = P*u_imag;

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
    
    [rhsu_real, rhsu_imag] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, u_real, ...
        u_imag, k, kb, gamma, sigmaTPA, sigma, P, W, Px, Py, Sx, Sy, Wx, Wy, xloc, yloc);
    u_real1 = u_real + c11*dt*rhsu_real;    
    u_imag1 = u_imag + c11*dt*rhsu_imag;
    
    [rhsu_real, rhsu_imag] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, u_real1, ...
        u_imag1, k, kb, gamma, sigmaTPA, sigma, P, W, Px, Py, Sx, Sy, Wx, Wy, xloc, yloc);
    u_real2 = c20*u_real + c21*u_real1 + c22*dt*rhsu_real;  
    u_imag2 = c20*u_imag + c21*u_imag1 + c22*dt*rhsu_imag;
    
    
    [rhsu_real, rhsu_imag] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, u_real2, ...
        u_imag2, k, kb, gamma, sigmaTPA, sigma, P, W, Px, Py, Sx, Sy, Wx, Wy, xloc, yloc);
    u_real3 = c30*u_real + c32*u_real2 + c33*dt*rhsu_real;   
    u_imag3 = c30*u_imag + c32*u_imag2 + c33*dt*rhsu_imag;
    
    
    [rhsu_real3, rhsu_imag3] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, u_real3, ...
        u_imag3, k, kb, gamma, sigmaTPA, sigma, P, W, Px, Py, Sx, Sy, Wx, Wy, xloc, yloc);
    u_real4 = c40*u_real + c43*u_real3 + c44*dt*rhsu_real3;   
    u_imag4 = c40*u_imag + c43*u_imag3 + c44*dt*rhsu_imag3;
    
    
    [rhsu_real, rhsu_imag]   = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, u_real4, ...
        u_imag4, k, kb, gamma, sigmaTPA, sigma, P, W, Px, Py, Sx, Sy, Wx, Wy, xloc, yloc);
    u_real  = c52*u_real2 + c53*u_real3 + c53_1*dt*rhsu_real3 + c54*u_real4 + c55*dt*rhsu_real; 
    u_imag  = c52*u_imag2 + c53*u_imag3 + c53_1*dt*rhsu_imag3 + c54*u_imag4 + c55*dt*rhsu_imag;
    
    ureal_ret(:,:,it) = P*u_real;
    uimag_ret(:,:,it) = P*u_imag;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error check

err_ureal = 0;
err_uimag = 0;
ucloc = zeros(((qx+1)*(qy+1)),(Nx*Ny));
ucloc_real = ucloc;
ucloc_imag = ucloc;
utrue_real = ucloc;
utrue_imag = ucloc;

for i=1:Nx
    for j=1:Ny
        [uloc] = solnew(xloc(:,i),yloc(:,j),T);
        
        utrue_real(:,(i-1)*Ny+j) = (uloc(1,:)');
        utrue_imag(:,(i-1)*Ny+j) = (uloc(2,:)');
        
        ucloc_real(:,(i-1)*Ny+j) = (P*u_real(:,(i-1)*Ny+j));
        ucloc_imag(:,(i-1)*Ny+j) = (P*u_imag(:,(i-1)*Ny+j));
        
        err_ureal = err_ureal + (hx(i)/2)*(hy(j)/2)*((ucloc_real(:,(i-1)*Ny+j) - ...
            utrue_real(:,(i-1)*Ny+j)))'*(W)*(ucloc_real(:,(i-1)*Ny+j) - utrue_real(:,(i-1)*Ny+j));

        err_uimag = err_uimag + (hx(i)/2)*(hy(j)/2)*((ucloc_imag(:,(i-1)*Ny+j) - ...
            utrue_imag(:,(i-1)*Ny+j)))'*(W)*(ucloc_imag(:,(i-1)*Ny+j) - utrue_imag(:,(i-1)*Ny+j));
    end
end

fprintf('u real L2 error = %4.3e \n',sqrt(err_ureal));
fprintf('u imag L2 error = %4.3e \n',sqrt(err_uimag));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot the projection

plot_ureal = (P*u_real);
plot_ureal_reshape = zeros((qx+1)*(qy+1),Nx*Ny);

plot_uimag = (P*u_imag);
plot_uimag_reshape = zeros((qx+1)*(qy+1),Nx*Ny);

plot_intensity = (P*u_real).^2 + (P*u_imag).^2;
plot_intensity_reshape = zeros((qx+1)*(qy+1),Nx*Ny);

for j=1:Ny
    for i=1:Nx
        plot_ureal_reshape(:,(j-1)*Ny+i) = plot_ureal(:,(i-1)*Ny+j);
        plot_uimag_reshape(:,(j-1)*Ny+i) = plot_uimag(:,(i-1)*Ny+j);
        plot_intensity_reshape(:,(j-1)*Ny+i) = plot_intensity(:,(i-1)*Ny+j);
    end
end

plot_z = zeros((qx+1)*Nx,(qy+1)*Ny);
plot_zimag = zeros((qx+1)*Nx,(qy+1)*Ny);
plot_zintensity = zeros((qx+1)*Nx,(qy+1)*Ny);

for j=1:Ny
    for i=1:Nx
        for k=1:(qy+1)
            new = plot_ureal_reshape(k,(j-1)*Nx+i);
            newimag = plot_uimag_reshape(k,(j-1)*Nx+i);
            newintensity = plot_intensity_reshape(k,(j-1)*Nx+i);
            
            for l=2:(qx+1)
                new(l) = plot_ureal_reshape(k+(l-1)*(qy+1),(j-1)*Nx+i);
                newimag(l) = plot_uimag_reshape(k+(l-1)*(qy+1),(j-1)*Nx+i);
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


function [u] = solnew(xloc,yloc,t)

dim_x = size(xloc);
len_x = dim_x(1); % qx+1
dim_y = size(yloc);
len_y = dim_y(1); % qy+1

u = zeros(2,len_x*len_y);

for i=1:len_x
    u(1,((i-1)*len_y+1):(i*len_y)) = cos(2*pi*(xloc(i)+yloc+t))';
    u(2,((i-1)*len_y+1):(i*len_y)) = sin(2*pi*(xloc(i)+yloc+t))';
end

end


function [ut_real, ut_imag] = schrodinger_update(qx, qy, hx, hy, Nx, Ny, Vb_x, Vb_y, Uxb_x, Uxb_y, u_real, ...
    u_imag, k, kb, gamma, sigmaTPA, sigma, P, W, Px, Py, Sx, Sy, Wx, Wy, xloc, yloc)

ut_real = zeros((qx+1)*(qy+1),(Nx*Ny));
ut_imag = zeros((qx+1)*(qy+1),(Nx*Ny));

for i=1:Nx
    for j=1:Ny
        
        % local value coressponds to (i,j)
        ureal_local = u_real(:, (i-1)*Ny + j);
        uimag_local = u_imag(:, (i-1)*Ny + j);
        
        k_local = k(:, (i-1)*Ny + j);
        gamma_local = gamma(:, (i-1)*Ny + j);
        sigmaTPA_local = sigmaTPA(:, (i-1)*Ny + j);
        sigma_local = sigma(:, (i-1)*Ny + j);
        
        % assemble ut
        Mu1 = (kron(Sx,Py)')*diag(1./(2.*k_local))*W*(kron(Sx,Py));
        Mu1 = Mu1 + (kron(Px,Sy)')*diag(1./(2.*k_local))*W*(kron(Px,Sy));
        
        fu2 = (P*ureal_local).^2 + (P*uimag_local).^2;
        Mu2 = (hx(1)/2)*(hy(1)/2)*(P')*diag(fu2)*diag(gamma_local)*W*P;
        
        Mu3 = (1/2)*(hx(1)/2)*(hy(1)/2)*(P')*diag(fu2)*diag(sigmaTPA_local)*W*P;
        
        Mu4 = (1/2)*(hx(1)/2)*(hy(1)/2)*(P')*diag(sigma_local)*W*P;
        
        ut_real(:,(i-1)*Ny+j) = (2/hx(1))*(2/hy(1))*( Mu1*uimag_local - Mu2*uimag_local ...
                               - Mu3*ureal_local - Mu4*ureal_local);
        ut_imag(:,(i-1)*Ny+j) = (2/hx(1))*(2/hy(1))*(-Mu1*ureal_local + Mu2*ureal_local ...
                               - Mu3*uimag_local - Mu4*uimag_local);
        
        % central flux: star = (local+out)/2
        % periodic BC
        % (i-1)*Ny + j
        if (j==1) % out corresponds to j == Ny (i-1)*Ny + Ny south side
           ureal_out = u_real(:, (i-1)*Ny + Ny);
           uimag_out = u_imag(:, (i-1)*Ny + Ny);
           kb_s = kb(1:(qx+1), (i-1)*Ny + Ny);
        else % out corresponds to j == j-1 (i-1)*Ny + j - 1
           ureal_out = u_real(:, (i-1)*Ny + j-1);
           uimag_out = u_imag(:, (i-1)*Ny + j-1);
           kb_s = kb(1:(qx+1), (i-1)*Ny + j-1);
        end
        
        Fu_real =         - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Vb_y(1,:)')*Uxb_y(1,:)) * (-uimag_local);
        Fu_real = Fu_real - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Vb_y(1,:)')*Uxb_y(2,:)) * (-uimag_out);
        Fu_real = Fu_real - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Uxb_y(1,:)')*Vb_y(1,:)) * (-uimag_local);
        Fu_real = Fu_real - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Uxb_y(1,:)')*Vb_y(2,:)) * ( uimag_out);
        
        Fu_imag =         - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Vb_y(1,:)')*Uxb_y(1,:)) * ( ureal_local);
        Fu_imag = Fu_imag - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Vb_y(1,:)')*Uxb_y(2,:)) * ( ureal_out);
        Fu_imag = Fu_imag - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Uxb_y(1,:)')*Vb_y(1,:)) * ( ureal_local);
        Fu_imag = Fu_imag - (1/2) * kron((Px')*diag(1./kb_s)*diag(Wx)*Px, (Uxb_y(1,:)')*Vb_y(2,:)) * (-ureal_out);
        
        if (j==Ny) % out corresponds to j == 1 north side
           ureal_out = u_real(:, (i-1)*Ny + 1);
           uimag_out = u_imag(:, (i-1)*Ny + 1);
           kb_n = kb(1:(qx+1),(i-1)*Ny + 1);
        else % out corresponds to j + 1
           ureal_out = u_real(:, (i-1)*Ny + j+1);
           uimag_out = u_imag(:, (i-1)*Ny + j+1);
           kb_n = kb(1:(qx+1),(i-1)*Ny + j+1);
        end
        
        Fu_real = Fu_real + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Vb_y(2,:)')*Uxb_y(2,:)) * (-uimag_local);
        Fu_real = Fu_real + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Vb_y(2,:)')*Uxb_y(1,:)) * (-uimag_out);
        Fu_real = Fu_real + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Uxb_y(2,:)')*Vb_y(2,:)) * (-uimag_local);
        Fu_real = Fu_real + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Uxb_y(2,:)')*Vb_y(1,:)) * ( uimag_out);
        
        Fu_imag = Fu_imag + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Vb_y(2,:)')*Uxb_y(2,:)) * ( ureal_local);
        Fu_imag = Fu_imag + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Vb_y(2,:)')*Uxb_y(1,:)) * ( ureal_out);
        Fu_imag = Fu_imag + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Uxb_y(2,:)')*Vb_y(2,:)) * ( ureal_local);
        Fu_imag = Fu_imag + (1/2) * kron((Px')*diag(1./kb_n)*diag(Wx)*Px, (Uxb_y(2,:)')*Vb_y(1,:)) * (-ureal_out);
        
        if (i==1) % out corresponds to i = Nx west side
           ureal_out = u_real(:, (Nx-1)*Ny + j);
           uimag_out = u_imag(:, (Nx-1)*Ny + j);
           kb_w = kb((qx+2):(qx+1+qy+1),(Nx-1)*Ny + j);
        else % out corresponds to i - 1
           ureal_out = u_real(:, (i-2)*Ny + j);
           uimag_out = u_imag(:, (i-2)*Ny + j);
           kb_w = kb((qx+2):(qx+1+qy+1),(i-2)*Ny + j);
        end
        
        Fu_real = Fu_real - (1/2) * kron((Vb_x(1,:)')*Uxb_x(1,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * (-uimag_local);
        Fu_real = Fu_real - (1/2) * kron((Vb_x(1,:)')*Uxb_x(2,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * (-uimag_out);
        Fu_real = Fu_real - (1/2) * kron((Uxb_x(1,:)')*Vb_x(1,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * (-uimag_local);
        Fu_real = Fu_real - (1/2) * kron((Uxb_x(1,:)')*Vb_x(2,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * ( uimag_out);
        
        Fu_imag = Fu_imag - (1/2) * kron((Vb_x(1,:)')*Uxb_x(1,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * ( ureal_local);
        Fu_imag = Fu_imag - (1/2) * kron((Vb_x(1,:)')*Uxb_x(2,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * ( ureal_out);
        Fu_imag = Fu_imag - (1/2) * kron((Uxb_x(1,:)')*Vb_x(1,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * ( ureal_local);
        Fu_imag = Fu_imag - (1/2) * kron((Uxb_x(1,:)')*Vb_x(2,:), (Py')*diag(1./kb_w)*diag(Wy)*Py) * (-ureal_out);
        
        if (i==Nx) % out corresponds to i = 1 east side
           ureal_out = u_real(:, j);
           uimag_out = u_imag(:, j);
           kb_e = kb((qx+2):(qx+1+qy+1),j);
        else % out corresponds to i + 1
           ureal_out = u_real(:, i*Ny + j);
           uimag_out = u_imag(:, i*Ny + j);
           kb_e = kb((qx+2):(qx+1+qy+1),i*Ny + j);
        end
        
        Fu_real = Fu_real + (1/2) * kron((Vb_x(2,:)')*Uxb_x(2,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * (-uimag_local);
        Fu_real = Fu_real + (1/2) * kron((Vb_x(2,:)')*Uxb_x(1,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * (-uimag_out);
        Fu_real = Fu_real + (1/2) * kron((Uxb_x(2,:)')*Vb_x(2,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * (-uimag_local);
        Fu_real = Fu_real + (1/2) * kron((Uxb_x(2,:)')*Vb_x(1,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * ( uimag_out);
        
        Fu_imag = Fu_imag + (1/2) * kron((Vb_x(2,:)')*Uxb_x(2,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * ( ureal_local);
        Fu_imag = Fu_imag + (1/2) * kron((Vb_x(2,:)')*Uxb_x(1,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * ( ureal_out);
        Fu_imag = Fu_imag + (1/2) * kron((Uxb_x(2,:)')*Vb_x(2,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * ( ureal_local);
        Fu_imag = Fu_imag + (1/2) * kron((Uxb_x(2,:)')*Vb_x(1,:), (Py')*diag(1./kb_e)*diag(Wy)*Py) * (-ureal_out);
        
        Fu_real = Fu_real/2;
        Fu_imag = Fu_imag/2;
        
        ut_real(:,(i-1)*Ny+j) = ut_real(:,(i-1)*Ny+j) + (2/hx(i))*(2/hy(j))*Fu_real;
        ut_imag(:,(i-1)*Ny+j) = ut_imag(:,(i-1)*Ny+j) + (2/hx(i))*(2/hy(j))*Fu_imag;
        
    end
end

end

