using FourierFlows, Printf, Random, FFTW

using Random: seed!

import GeophysicalFlows.TwoDNavierStokes
import GeophysicalFlows.TwoDNavierStokes: energy, enstrophy
import GeophysicalFlows: peakedisotropicspectrum

dev = CPU()     # Device (CPU/GPU)

# Parameters
  n = 256
  L = 2π
 nν = 2
  ν = 0.0
 dt = 1e-3
nstepsjit = 10 # call stepforward for nstepsjit to force compilation
nsteps = 5000
 nothingfunction(args...) = nothing


gr = TwoDGrid(dev, n, L, n, L; nthreads=1)
pr = TwoDNavierStokes.Params(ν, nν)
vs = TwoDNavierStokes.Vars(dev, gr)
eq = TwoDNavierStokes.Equation(pr, gr)

prob = FourierFlows.Problem(eq, "FilteredAB3", dt, gr, vs, pr, dev)
filter = FourierFlows.makefilter(prob.eqn)

# some aliases
sol, cl, vs, gr = prob.sol, prob.clock, prob.vars, prob.grid

x, y = gridpoints(gr)

# Initial condition closely following pyqg barotropic example
# that reproduces the results of the paper by McWilliams (1984)
seed!(1234)
k0, E0 = 6, 0.5
zetai  = peakedisotropicspectrum(gr, k0, E0, mask=filter)
TwoDNavierStokes.set_zeta!(prob, zetai)

for j=1:nstepsjit  #just-in-time compilation
  stepforward!(prob)
end

TwoDNavierStokes.set_zeta!(prob, zetai)

startwalltime = time()
while cl.step < nsteps
  stepforward!(prob)
end

println(round((time()-startwalltime)/(nsteps-1)*1000, digits=3), " ms per time-step")

# using PyPlot
# TwoDNavierStokes.updatevars!(prob)
# figure(figsize=(10, 4))
# subplot(121)
# pcolormesh(x, y, zetai)
# xlabel("x")
# ylabel("y")
# title("vorticity @ t=0")
# axis("square")
# subplot(122)
# pcolormesh(x, y, vs.zeta)
# axis("square")
# xlabel("x")
# ylabel("y")
# title("vorticity @ t="*string(round(cl.t, digits=2)))
# savefig("GeophysicalFlows_n"*string(n)*".png", dpi=400)
