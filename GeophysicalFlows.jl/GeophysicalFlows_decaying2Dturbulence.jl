using FourierFlows, Printf

import GeophysicalFlows.TwoDNavierStokes

using Random: seed!
using Statistics: mean

dev = CPU()     # Device (CPU/GPU)

# Parameters
 n = 256
 L = 2π
 ν = 1e-4
nν = 1
dt = 1e-3
nsteps_jit = 10 # call stepforward for nsteps_jit to force compilation
nsteps = 5000
stepper = "AB3"
floattype = Float64

grid2D = TwoDGrid(dev, n, L, n, L; nthreads=1, T=floattype)
params = TwoDNavierStokes.Params(ν, nν)
vars = TwoDNavierStokes.Vars(dev, grid2D)
equation = TwoDNavierStokes.Equation(params, grid2D)

prob = FourierFlows.Problem(equation, stepper, dt, grid2D, vars, params, dev)

# Random initial condition
seed!(1234)
zeta_initial = randn(floattype, (n, n))
zeta_initial = zeta_initial .- mean(zeta_initial) # make sure initial condition has zero mean
TwoDNavierStokes.set_zeta!(prob, zeta_initial)

for j=1:nsteps_jit  #just-in-time compilation
  stepforward!(prob)
end

TwoDNavierStokes.set_zeta!(prob, zeta_initial)

startwalltime = time()
while prob.clock.step < nsteps
  stepforward!(prob)
  dealias!(prob.sol, prob.grid)
end

println(round((time()-startwalltime)/(nsteps-1)*1000, digits=3), " ms per time-step")

# using PyPlot
# x, y = gridpoints(prob.grid)
# TwoDNavierStokes.updatevars!(prob)
# figure(figsize=(10, 4))
# subplot(121)
# pcolormesh(x, y, zeta_initial)
# xlabel("x")
# ylabel("y")
# title("vorticity @ t=0")
# axis("square")
# subplot(122)
# pcolormesh(x, y, prob.vars.zeta)
# axis("square")
# xlabel("x")
# ylabel("y")
# title("vorticity @ t="*string(round(prob.clock.t, digits=2)))
# savefig("GeophysicalFlows_n"*string(n)*".png", dpi=400)
