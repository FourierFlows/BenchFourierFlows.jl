clear all;

global invKsq K L Ksq nu nnu filter

maxNumCompThreads(1);
GPU = 0; % set =1 to run on GPU

nx = 256;
ny = nx;

Lx = 2*pi;
Ly = Lx;

nu = 0.0;
nnu = 2;

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

% creat a filter in wavenumber space;
% a is chosen so that the energy at the largest nondim
% wavenumber K*dx be zero whithin machine double precision
s=4;
Kmax = ny/2; Kmax_s=Kmax*dy;
kcut = 2/3*Kmax;kcut_s=kcut*dy;
a = -log(1e-15)/(Kmax_s-kcut_s)^s * dy^s;
Kmag=sqrt(Ksq);

filter = ones(ny, nx).*abs(Kmag<=ny/3) + exp(-a*(Kmag-kcut).^s).*abs(Kmag>ny/3);

psih0 = (Ksq.*(1 + (Kmag/6).^4)).^(-1/2) .* (randn(ny, nx) + 1i*randn(ny, nx));
psih0(1, 1)=0;
psih0 = fft2(real(ifft2(psih0)));
Etemp = sum(0.5*Ksq(:).*abs(psih0(:)).^2/(nx*ny)^2);
psih0 = psih0*sqrt(0.5/Etemp);
qh0 = -Ksq.*psih0.*filter;
qh = qh0;

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
% pcolor2(X, Y, real(ifft2(qh0)))
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
% print('matlab_n256.png', '-dpng', '-r400'); 