import numpy as np
import time, logging

from mpi4py import MPI
from dedalus import public as de

from numpy import pi

logger = logging.getLogger(__name__)
minute = 60.0
hour = 60*minute

Lx, Ly = 2*pi, 2*pi
nx, ny = 256, 256
nu = 0.

xbasis = de.Fourier('x', nx, interval=(-Lx/2, Lx/2), dealias=1)
ybasis = de.Fourier('y', ny, interval=(-Ly/2, Ly/2), dealias=1)
domain = de.Domain([xbasis, ybasis], grid_dtype=np.float64)

x, y = domain.grid(0), domain.grid(1)

variables = ['q', 'psi']
problem = de.IVP(domain, variables=variables, time='t')

problem.parameters['nu'] = nu

problem.substitutions['J(a, b)'] = "dx(a)*dy(b) - dy(a)*dx(b)"
problem.substitutions['lap(a)'] = "d(a, x=2) + d(a, y=2)"

# Equations
problem.add_equation("dt(q) + nu*lap(lap(q)) = - J(psi, q)")
problem.add_equation("q - lap(psi) = 0", condition="(nx != 0) or (ny != 0)")
problem.add_equation("psi = 0", condition="(nx == 0) and (ny == 0)")


start_build_time = time.time()
solver = problem.build_solver(de.timesteppers.SBDF2)
logger.info('Solver built. (t = %f) ' %(time.time()-start_build_time))

q = solver.state['q']
psi = solver.state['psi']

def constructfilter(xbasis, ybasis):
    kx, ky = np.meshgrid(xbasis.wavenumbers, ybasis.wavenumbers, indexing='ij')
    nkx, nky = kx.shape
    sq = kx**2+ky**2
    cphi=.65*np.pi;
    filterfac=23.6;
    wv = np.sqrt((kx*Lx/nx)**2+(ky*Ly/ny)**2)
    filter = np.exp(-filterfac*(wv-cphi)**4)
    filter[wv<=cphi] = 1.
    return filter
    
filter = constructfilter(xbasis, ybasis)

# Initialize state variables
def peakedisotropicspectrum(xbasis, ybasis, domain, k0=6, energy0=0.5, seed=1234):
    kx, ky = np.meshgrid(xbasis.wavenumbers, ybasis.wavenumbers, indexing='ij')
    nkx, nky = kx.shape
    x, y = np.meshgrid(xbasis.grid(), ybasis.grid(), indexing='ij')
    nx, ny = x.shape
    
    modk = np.sqrt(kx**2+ky**2)
    psik = (modk**2 * (1 + (modk/k0)**4)+1e-14)**(-0.5)
    psik[modk==0] = 0.
    
    rand = np.random.RandomState(seed=seed)
    cshape = domain.dist.coeff_layout.global_shape(scales=1)
    slices = domain.dist.coeff_layout.slices(scales=1)
    
    phases = rand.standard_normal(cshape)[slices] + 1j*rand.standard_normal(cshape)[slices]
    psih = phases*psik
    Ein = (modk**2 * np.abs(psih)**2/(nx*ny)**2).sum()
    psih = psih*np.sqrt(energy0/Ein)
    q, psi = np.fft.irfft2(-modk**2*psih, (nx, ny)), np.fft.irfft2(psih, (nx, ny))
    return q, psi

q['g'], psi['g'] = peakedisotropicspectrum(xbasis, ybasis, domain, k0=6, energy0=0.5)
q['c'], psi['c'] = q['c']*filter, psi['c']*filter

# Initial timestep
dt = 0.01

# Integration parameters
solver.stop_sim_time = np.inf
solver.stop_wall_time = np.inf
solver.stop_iteration = 500
log_cadence=500

def time_to_log(log_cadence):
    (solver.iteration-1) % log_cadence == 0

# Main loop
try:
    logger.info('Starting loop')
    start_run_time = time.time()
    while solver.ok:
        solver.step(dt)
        q['c'] = q['c']*filter
        psi['c'] = psi['c']*filter
        if time_to_log(log_cadence): 
            log(logger, dt)
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Time per time-step:  %.3f ms' %((end_run_time-start_run_time)/solver.iteration*1000))
    # logger.info(
    #     'Run time: %f cpu-hr' %((end_run_time-start_run_time)/hour * domain.dist.comm_cart.size))


