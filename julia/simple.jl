using LinearAlgebra: mul!, ldiv!
# using PyPlot
using FFTW
# using CuArrays

ArrayType=Array # or CuArrays
T = Float32
n = 256
dt = T(1e-3)
nsteps = 5000
 
 nx, ny = n, n
 Lx, Ly = 2π, 2π
 nk, nl = nx, ny
nkr = Int(nx/2+1)
nu, nnu = 0.0, 2

x = reshape(ArrayType(range(-Lx/2, step=Lx/nx, length=nx)), (nx, 1))
y = reshape(ArrayType(range(-Ly/2, step=Ly/ny, length=ny)), (1, ny))

X, Y = x .+ 0*y, 0*x .+ y

 # Wavenubmer grid
i₁ = 0:Int(nx/2)
j₁ = 0:Int(ny/2)
j₂ = Int(-ny/2+1):-1

kr = ArrayType(reshape(2π/Lx*cat(i₁, dims=1), (nkr, 1)))
 l = ArrayType(reshape(2π/Ly*cat(j₁, j₂, dims=1), (1, nl)))

   Krsq = @. kr^2 + l^2
invKrsq = @. 1/Krsq
invKrsq[1, 1] = 0

nthreads=1

# FFT plans
FFTW.set_num_threads(nthreads)
rfftplan = plan_rfft(ArrayType{T, 2}(undef, nx, ny))

function zeroarray(dims; T=Float64, A=ArrayType)
  a = A{T}(undef, dims...)
  a .= 0
  return a
end

u = zeroarray((nx, ny), T=T)
v = zeroarray((nx, ny), T=T)
q = zeroarray((nx, ny), T=T)

uh = zeroarray((nkr, nl), T=Complex{T})
vh = zeroarray((nkr, nl), T=Complex{T})
qh = zeroarray((nkr, nl), T=Complex{T})

RHS   = zeroarray((nkr, nl), T=Complex{T})
RHS₋₁ = zeroarray((nkr, nl), T=Complex{T})
RHS₋₂ = zeroarray((nkr, nl), T=Complex{T})


function makefilter(kr, l, dx, dy, Tsol; A=ArrayType, order=4, innerK=0.65, outerK=1)
  K = @. (sqrt((kr*dx/π)^2 + (l*dy/π)^2))
  K = Array(K)
  decay = 15*log(10) / (outerK-innerK)^order # decay rate for filtering function
  filt = @. exp( -decay*(K-innerK)^order )
  filt[real.(K) .< innerK] .= 1
  Tsol(filt)
end

filterq = makefilter(kr, l, Lx/nx, Ly/ny, typeof(qh))

function peakedisotropicspectrum(kpeak::Real, E0::Real; T=Float64, A=ArrayType, mask=ones(size(Krsq)))
  k0 = kpeak*2π/Lx
  modk = sqrt.(Krsq)
  ψk = A(zeros(T, (nk, nl)))
  ψk = @. (modk^2 * (1 + (modk/k0)^4))^(-0.5)
  ψk[1, 1] = 0.0
  phases = randn(Complex{T}, size(Krsq))
  phases_real, phases_imag = real.(phases), imag.(phases)
  phases = A(phases_real) + im*A(phases_imag)
  ψh = @. phases*ψk
  ψh = ψh.*A(mask)
  Ein = real(sum(Krsq.*abs2.(ψh)/(nx*ny)^2))
  ψh = ψh*sqrt(E0/Ein)
  q = A(-irfft(Krsq.*ψh, nx))
end

kpeak, E0 = 6, 0.5
q = typeof(q)(peakedisotropicspectrum(kpeak, E0))
qh = typeof(qh)(rfft(q))

function RHS_eq!(RHS, qh)
  uh = @.  im * l  * invKrsq * qh
  vh = @. -im * kr * invKrsq * qh

  ldiv!(u, rfftplan, uh)
  ldiv!(v, rfftplan, vh)
  ldiv!(q, rfftplan, qh)
  
  @. u *= q
  @. v *= q
  
  mul!(uh, rfftplan, u)
  mul!(vh, rfftplan, v)
  
  @. RHS = - im*kr*uh - im*l*vh - nu*Krsq^nnu*qh
end

function ForwardEulerTimeStep!(qh, dt, filterq)
  RHS_eq!(RHS, qh)
  @. qh += RHS*dt
  @. qh *= filterq
end

function AB3TimeStep!(qh, dt, RHS₋₁, RHS₋₂, filterq)
  RHS_eq!(RHS, qh)
  @. qh += (23*RHS - 16*RHS₋₁ + 5*RHS₋₂)*dt/12
  @. qh *= filterq
  @. RHS₋₂ = RHS₋₁
  @. RHS₋₁ = RHS
end

for i=1:3
  # call this function to force JIT compilation
  AB3TimeStep!(qh, dt, RHS₋₁, RHS₋₂, filterq)
end
qh = typeof(qh)(rfft(q))

for i = 1:2
  ForwardEulerTimeStep!(qh, dt, filterq)
end

startwalltime = time()
for i = 3:nsteps
  AB3TimeStep!(qh, dt, RHS₋₁, RHS₋₂, filterq)
end
elapsed = (time()-startwalltime)/(nsteps-2)*1000

println(round(elapsed, digits=3), " ms per time-step")
