using FourierFlows, Printf, Random, FFTW

using Random: seed!

import GeophysicalFlows.TwoDTurb
import GeophysicalFlows.TwoDTurb: energy, enstrophy
import GeophysicalFlows: peakedisotropicspectrum

dev = CPU()     # Device (CPU/GPU)

# Parameters
  n = 512
  L = 2π
 nν = 2
  ν = 0.0
 dt = 1e-2
nsteps = 500
 nsubs = 500
 T = Float64
 nothingfunction(args...) = nothing


gr = TwoDGrid(dev, n, L, n, L; T=T, nthreads=4, effort=FFTW.PATIENT)
pr = TwoDTurb.Params{T}(ν, nν, 0, 0, nothingfunction)
vs = TwoDTurb.Vars(dev, gr)
eq = TwoDTurb.Equation(pr, gr)
prob = FourierFlows.Problem(eq, "FilteredAB3", dt, gr, vs, pr, dev)

sol, cl, vs, gr, filter = prob.sol, prob.clock, prob.vars, prob.grid, prob.timestepper.filter
prob = FourierFlows.Problem(eq, "AB3", dt, gr, vs, pr, dev)
sol, cl, vs, gr = prob.sol, prob.clock, prob.vars, prob.grid
x, y = gridpoints(gr)

# Initial condition closely following pyqg barotropic example
# that reproduces the results of the paper by McWilliams (1984)
seed!(1234)
k0, E0 = 6, 0.5
zetai  = peakedisotropicspectrum(gr, k0, E0, mask=filter)
TwoDTurb.set_zeta!(prob, zetai)

# # Create Diagnostic -- energy and enstrophy are functions imported at the top.
# E = Diagnostic(energy, prob; nsteps=nsteps)
# Z = Diagnostic(enstrophy, prob; nsteps=nsteps)
# diags = [E, Z] # A list of Diagnostics types passed to "stepforward!" will
# # be updated every timestep.

startwalltime = time()
while cl.step < nsteps
  stepforward!(prob, nsubs)
  
  # Message
  log = @sprintf("step: %04d, t: %d, τ: %.2f min",
    cl.step, cl.t, (time()-startwalltime)/60)

  println(log)
end
println((time()-startwalltime))
println("finished")
