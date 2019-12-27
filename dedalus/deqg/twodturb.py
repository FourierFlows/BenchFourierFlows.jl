import numpy as np
import time, logging

from mpi4py import MPI
from dedalus import public as de
from numpy import pi

from .quasigeostrophy import QGModel
from .utils import add_parameters, add_first_derivative_substitutions, bind_parameters
from .utils import bind_state_variables

logger = logging.getLogger(__name__)

class TwoDTurbModel(QGModel):
    """
    A model for 2D turbulence.

    Args
    ----
        nx : (int)
            Grid resolution in :math:`x`

        ny : (int)
            Grid resolution in :math:`y`

        Lx : (float)
            Domain extent in :math:`x`

        Ly : (float)
            Domain extent in :math:`y`

        ν : (float)
            Hyperviscosity

        **params : (any)
            Additional parameters to be added to the dedalus problem.
    """
    def __init__(self,
        nx = 170,
        ny = 170,
        Lx = 2*pi,
        Ly = 2*pi,
        ν = 0.0,
        **params
        ):

        # Create bases and domain
        self.xbasis = xbasis = de.Fourier('x', nx, interval=(-Lx/2, Lx/2), dealias=3/2)
        self.ybasis = ybasis = de.Fourier('y', ny, interval=(-Ly/2, Ly/2), dealias=3/2)
        self.domain = domain = de.Domain([xbasis, ybasis], grid_dtype=np.float64)

        self.variables = variables = ['q', 'ψ']
        self.problem = problem = de.IVP(domain, variables=variables, time='t')

        add_parameters(problem, ν=ν, **params)
        bind_parameters(self, ν=ν, **params)

        add_first_derivative_substitutions(problem, ['q', 'ψ'], ['x', 'y'])

        problem.substitutions['J(a, b)'] = "dx(a)*dy(b) - dy(a)*dx(b)"
        problem.substitutions['lap(a)'] = "dx(dx(a)) + dy(dy(a))"
        
        # Equations
        problem.add_equation("dt(q) + ν*lap(lap(q)) = - J(ψ, q)")
        problem.add_equation("q - lap(ψ) = 0", condition="(nx != 0) or (ny != 0)")
        problem.add_equation("ψ = 0", condition="(nx == 0) and (ny == 0)")
        
        self.x = domain.grid(0)
        self.y = domain.grid(1)

    def build_solver(self, timestepper='RK443'):
        """Build a dedalus solver for the model with `timestepper`.

        Args
        ----
            timestepper : The name of the timestepper to be used by the solver.
        """
        QGModel.build_solver(self, timestepper=timestepper)
        bind_state_variables(self, self.variables)
