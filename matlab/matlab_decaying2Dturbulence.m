clear all;

global invKsq K L Ksq nu nnu

maxNumCompThreads(1);
GPU = 0; % set equal to 1 to run on GPU

nx = 256;
ny = nx;

Lx = 2*pi;
Ly = Lx;

nu = 1e-4;
nnu = 1;

dt = 1e-3;
nsteps = 5000;

T = zeros(nsteps+1, 1);

dx = Lx/nx;
dy = Ly/ny;

x = -Lx/2:dx:Lx/2-dx;
y = -Ly/2:dy:Ly/2-dy;

k = 2*pi/Lx * [0:nx/2-1 -nx/2:-1];
l = 2*pi/Ly * [0:nx/2-1 -nx/2:-1];

[X, Y] = meshgrid(x, y);
[K, L] = meshgrid(k, l);

Ksq = K.^2 + L.^2;
invKsq = 1 ./ Ksq;
invKsq(Ksq==0)=0;

% Random initial condition
q0 = randn(ny, nx);
q0 = q0 - mean(q0);  % make sure initial condition has zero mean
qh = fft2(q0);

if GPU==1
    qh = gpuArray(qh); K = gpuArray(K); L = gpuArray(L); Ksq = gpuArray(Ksq);
    invKsq = gpuArray(invKsq); filter = gpuArray(filter);
end

% first two time-steps using Forward Euler
T(2) = dt;
RHSpp = RHS_eq(qh);
qh = ForwardEulertimestep(qh, dt);

T(3) = 2*dt;
RHSp = RHS_eq(qh);
qh = ForwardEulertimestep(qh, dt);

tic
for it=3:nsteps
    T(it+1) = it*dt;
    [qh, RHSp, RHSpp] = AB3timestep(qh, RHSp, RHSpp, dt);
end
elapsed = toc;

disp(['Time per time-step: ', num2str(elapsed/(nsteps-2)*1000, '%1.3f'), ' ms']);

if GPU==1
    qh = gather(qh); K = gather(K); L = gather(L); Ksq = gather(Ksq);
    invKsq = gather(invKsq); filter = gather(filter);
end

% figure()
% set(gcf,'units','inches','position',[3, 3, 10, 4])
% subplot(1,2,1)
% pcolor2(X, Y, q0)
% set(gca, 'fontsize', 12, 'Box', 'on', 'layer', 'top')
% xticks(-3:3); yticks(-3:3);
% xlabel('$x$', 'fontsize', 16, 'Interpreter', 'latex')
% ylabel('$y$', 'fontsize', 16, 'Interpreter', 'latex')
% title('vorticity @ $t=0$', 'fontsize', 16, 'Interpreter', 'latex')
% axis square
% subplot(1,2,2)
% pcolor2(X, Y, real(ifft2(qh)))
% set(gca, 'fontsize', 12, 'Box', 'on', 'layer', 'top')
% xticks(-3:3); yticks(-3:3);
% xlabel('$x$', 'fontsize', 16, 'Interpreter', 'latex')
% ylabel('$y$', 'fontsize', 16, 'Interpreter', 'latex')
% title(['vorticity @ $t=' num2str(T(end)) '$'], 'fontsize', 16, 'Interpreter', 'latex')
% axis square;
% set(gcf,'PaperPositionMode','auto')
% print(['matlab_n' num2str(nx) '.png'], '-dpng', '-r400');