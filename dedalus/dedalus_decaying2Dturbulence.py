import sys; sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np

from numpy import pi
from deqg import TwoDTurbModel, mpiprint, random_noise, peakedisotropicspectrum

mpiprint("executable: {}".format(sys.executable))

# Calculation of ν (4th order hyperviscosity):
#
# >> ν = C * qmax / k^4, k=(dissipation scale).
#
# with C=10, k=(nx/2)/Lx=16, qmax=1e-4,
# => ν = 1.5e-8.

Lx, Ly = 2*pi, 2*pi
nx, ny = 70, 70

model = TwoDTurbModel(
    nx = nx,
    ny = ny,
    Lx = Lx,
    Ly = Ly,
    ν = 0.e-10,
    pi = pi,
)

model.build_solver()

# initial condition...
k0, E0 = 6, 0.5
model.set_fields(
    q = peakedisotropicspectrum(model.domain, k0=6, energy0=0.5)
)

kx, ky = np.meshgrid(model.xbasis.wavenumbers, model.ybasis.wavenumbers, indexing='ij')
nkx, nky = kx.shape
sq = kx**2+ky**2
cphi=.65*np.pi;
filterfac=23.6;
wv = np.sqrt((kx*Lx/nx)**2+(ky*Ly/ny)**2)
filtr = np.exp(-filterfac*(wv-cphi)**4)
filtr[wv<=cphi] = 1.

X, Y = np.meshgrid(model.x, model.y)

model.q['c'] = model.q['c']*filtr
model.ψ['c'] = model.ψ['c']*filtr

dt = 0.01

fig, axs = plt.subplots(1, 2)
for i in range(200):
    model.stop_at(iteration=100)
    model.run(dt=dt, log_cadence=100)

    model.q.set_scales(1)
    model.ψ.set_scales(1)
    
    axs[0].clear()
    imq1 = axs[0].contourf(X, Y, model.q['g'].T)
    # fig.colorbar(imq1, cax=axs[0, 0])
    # axs[1, 0].cla()
    
    # fig.colorbar(imq2, cax=axs[1, 0])
    axs[1].cla()
    imψ1 = axs[1].contourf(X, Y, model.ψ['g'].T)
    # fig.colorbar(imψ1, cax=axs[0, 1])
    axs[1].contour(X, Y, model.ψ['g'].T, colors="k")
    # axs[1, 1].cla()
    plt.pause(0.1)

    # fac = 1e-5/np.max(np.abs(model.q['g']))
    # 
    # model.q['c'] = fac*model.q['c']
    # model.ψ['c'] = fac*model.ψ['c']
    # 
    mean_q = np.mean(model.q['g'])
    mpiprint("mean pv: {}".format(mean_q))
