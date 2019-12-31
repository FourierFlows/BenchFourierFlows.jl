using FourierFlows, Printf, Random, FFTW

using Random: seed!

import GeophysicalFlows.TwoDTurb
import GeophysicalFlows.TwoDTurb: energy, enstrophy
import GeophysicalFlows: peakedisotropicspectrum

dev = CPU()     # Device (CPU/GPU)

# Parameters
  n = 256
  L = 2π
 nν = 2
  ν = 0.0
 dt = 1e-2
nstepsjit = 5
nsteps = 500
 nsubs = 500
 T = Float64
 nothingfunction(args...) = nothing


gr = TwoDGrid(dev, n, L, n, L; T=T, nthreads=1, effort=FFTW.PATIENT)
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

for j=1:nstepsjit  #just-in-time compilation
  stepforward!(prob)
  dealias!(sol, gr)
end

startwalltime = time()
while cl.step < nsteps+nstepsjit
  stepforward!(prob)
  dealias!(sol, gr)
end
# # Message
# log = @sprintf("step: %04d, t: %d, τ: %.2f min",
#   cl.step, cl.t, (time()-startwalltime)/60)
# 
# println(log)
println(round((time()-startwalltime)/(nsteps-1)*1000, digits=3), " ms per time-step")
println("finished")
