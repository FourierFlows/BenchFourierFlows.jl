import numpy as np
import time, logging

from mpi4py import MPI
from dedalus import public as de

logger = logging.getLogger(__name__)

minute = 60.0
hour = 60*minute

class QGModel():
    def stop_at(self, sim_time=np.inf, wall_time=np.inf, iteration=np.inf):
        """Direct the solver to stop the simulation at the
        specified `sim_time, `wall_time`, or `iteration`.

        Args
        ----
            sim_time : float
                Simulation or model time at which the solver stops.

            wall_time : float
                Wall time (time recorded by a 'clock on the wall', external to the simulation)
                at which the solver stops.

            iteration : float
                Time iteration at which the solver stops.
        """
        self.solver.stop_sim_time = sim_time
        self.solver.stop_wall_time = wall_time
        self.solver.stop_iteration = self.solver.iteration + iteration

    def set_field(self, phi, gridvalue):
        """Set `phi` as `gridvalue`. Calculate derivatives of `phi`.

        Args
        ----
            phi : str
                The name of the field to be set.

            gridvalue : np.ndarray
                The grid values of phi.
        """

        field = getattr(self, phi)
        field['g'] = gridvalue

        for dim in ('x', 'y', 'z'):
            if hasattr(self, phi+dim):
                field.differentiate(dim, out=getattr(self, phi+dim))

    def set_fields(self, **fields):
        """Set the fields defined by `fields`."""
        for name, value in fields.items():
            self.set_field(name, value)

    def add_log_tasks(self, **tasks):
        """Add tasks to be logged during model run."""
        for name, task in tasks.items():
            try:
                self.log_tasks[name] = task
            except AttributeError:
                self.log_tasks = {name: task}

    def build_solver(self, timestepper='RK443'):
        """Build a dedalus solver for the model with `timestepper`.

        Args
        ----
            timestepper : The name of the timestepper to be used by the solver. 
        """
        detimestepper = getattr(de.timesteppers, timestepper)

        start_build_time = time.time()
        solver = self.problem.build_solver(detimestepper)
        logger.info('Solver built. (t = %f) ' %(time.time()-start_build_time))

        self.solver = solver

    def log(self, logger, dt):
        """Print informational messages about a model run."""
        logger.info('Iteration: %i, Time: %e, dt: %e' %(self.solver.iteration, self.solver.sim_time, dt))

        try:
            for name, task in self.log_tasks.items(): 
                logger.info("{}: {}".format(name, task(self)))
        except AttributeError:
            pass

    def time_to_log(self, log_cadence):
        """ Return True if it is a logging iteration."""
        return (self.solver.iteration-1) % log_cadence == 0

    def run(self, dt=1e-16, log_cadence=100, runlogger=logger):
        """Run a model.

        Args
        ----
            dt : float
                The time step.

            log_cadence : int
                How often the simulation logs output.

            runlogger : int
                The logger to use for logging output during the run.
        """

        try:
            runlogger.info('Starting loop')
            start_run_time = time.time()
            while self.solver.ok:
                self.solver.step(dt)
                if self.time_to_log(log_cadence): 
                    self.log(runlogger, dt)
        except:
            runlogger.error('Exception raised, triggering end of main loop.')
            raise
        finally:
            end_run_time = time.time()
            runlogger.info('Iterations: %i' %self.solver.iteration)
            runlogger.info('Sim end time: %f' %self.solver.sim_time)
            runlogger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
            runlogger.info(
                'Run time: %f cpu-hr' %((end_run_time-start_run_time) / hour * self.domain.dist.comm_cart.size))
