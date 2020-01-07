import numpy as np
import time, logging

from mpi4py import MPI
from dedalus import public as de
from dedalus.extras import plot_tools

from numpy import pi

# Parameters
Lx = Ly = 2 * pi
nx = ny = 256
dealias = 1
nu = 0.
dt = 0.01
stop_iteration = 500
startup_iterations = 10

logger = logging.getLogger(__name__)
minute = 60.0
hour = 60*minute

xbasis = de.Fourier('x', nx, interval=(-Lx/2, Lx/2), dealias=dealias)
ybasis = de.Fourier('y', ny, interval=(-Ly/2, Ly/2), dealias=dealias)
domain = de.Domain([xbasis, ybasis], grid_dtype=np.float64)

x, y = domain.grid(0), domain.grid(1)

variables = ['q', 'psi']
problem = de.IVP(domain, variables=variables, time='t')

problem.parameters['nu'] = nu

problem.substitutions['J(a, b)'] = "dx(a)*dy(b) - dy(a)*dx(b)"
problem.substitutions['lap(a)'] = "d(a, x=2) + d(a, y=2)"

# Equations
problem.add_equation("dt(q) + nu*lap(lap(q)) = - J(psi, q)", condition="(nx != 0) or (ny != 0)")
problem.add_equation("q - lap(psi) = 0")
problem.add_equation("psi = 0", condition="(nx == 0) and (ny == 0)")


start_build_time = time.time()
solver = problem.build_solver(de.timesteppers.SBDF2)
logger.info('Solver built. (t = %f) ' %(time.time()-start_build_time))

q = solver.state['q']
psi = solver.state['psi']

def constructfilter(domain):
    kx = domain.elements(0)
    ky = domain.elements(1)
    cphi = 0.65 * np.pi
    filterfac = 23.6
    wv = np.sqrt((kx*Lx/nx)**2+(ky*Ly/ny)**2)
    filter = np.exp(-filterfac*(wv-cphi)**4)
    filter[wv<=cphi] = 1.
    return filter

filter = constructfilter(domain)

# Initialize state variables
def peakedisotropicspectrum(domain, k0=6, energy0=0.5, seed=1234):
    # Wavenumbers
    kx = domain.elements(0)
    ky = domain.elements(1)
    modk = np.sqrt(kx**2+ky**2)
    # Isotropic spectrum
    psi = domain.new_field()
    psi['c'] = (modk**2 * (1 + (modk/k0)**4) + 1e-14)**(-0.5)
    psi['c'][modk == 0] = 0
    # Add random phases, globally initialized for parallel reproducibility
    rand = np.random.RandomState(seed=seed)
    cshape = domain.dist.coeff_layout.global_shape(scales=1)
    slices = domain.dist.coeff_layout.slices(scales=1)
    phases = rand.standard_normal(cshape)[slices] + 1j*rand.standard_normal(cshape)[slices]
    psi['c'] *= phases
    # Impose Hermitian symmetry
    psi['g']
    # Normalize energy
    u = psi.differentiate('x')
    v = psi.differentiate('y')
    Ein = (0.5 * (u*u + v*v)).evaluate().integrate()
    psi['g'] *= (energy0 / Ein['g'])**0.5
    return psi

psi['g'] = peakedisotropicspectrum(domain, k0=6, energy0=0.5)['g']
q['c'] = problem.namespace['lap'](psi).evaluate()['c']
psi['c'] *= filter
q['c'] *= filter
qi = q['g'].copy()

# Integration parameters
solver.stop_sim_time = np.inf
solver.stop_wall_time = np.inf
solver.stop_iteration = stop_iteration
log_cadence = stop_iteration

def time_to_log(log_cadence):
    (solver.iteration-1) % log_cadence == 0

# Main loop
try:
    logger.info('Starting loop')
    while solver.ok:
        if solver.iteration == startup_iterations:
            start_run_time = time.time()
        solver.step(dt)
        q['c'] = q['c']*filter
        # psi['c'] = psi['c']*filter
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
    logger.info('Time per time-step:  %.3f ms' %((end_run_time-start_run_time)/(solver.iteration - startup_iterations)*1000))
    # logger.info(
    #     'Run time: %f cpu-hr' %((end_run_time-start_run_time)/hour * domain.dist.comm_cart.size))

# Gather distributed snapshots
qf = q['g'].copy()
qi = domain.dist.comm.gather(qi, root=0)
qf = domain.dist.comm.gather(qf, root=0)

# Plot from root
if domain.dist.comm.rank == 0:
    import matplotlib.pyplot as plt
    qi = np.concatenate(qi, axis=1)
    qf = np.concatenate(qf, axis=1)
    X, Y = plot_tools.quad_mesh(xbasis.grid(1), ybasis.grid(1))
    plt.figure(figsize=(10, 4))
    plt.subplot(121)
    plt.pcolormesh(X, Y, qi.T)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("vorticity @ t=0")
    plt.axis("square")
    plt.subplot(122)
    plt.pcolormesh(X, Y, qf.T)
    plt.axis("square")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("vorticity @ t= %.2f" %solver.sim_time)
    plt.savefig("dedalus_n256.png", dpi=400)
