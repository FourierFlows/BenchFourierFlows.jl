# BenchFourierFlows.jl
This project benchmarks FourierFlows performance against other similar packages

Results for `nx=ny=512` grid-points:

On *CPU* (2.9 GHz Intel Core i7 Macbook Pro with 16 GB 2133 MHz LPDDR3, using a single thread)
- [GeophysicalFlows.jl v0.3.0](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/v0.3.0): **7.64 ms** per time-step
- [pyqg v0.2.0](https://github.com/pyqg/pyqg/tree/v0.2.0): **12.50 ms** per time-step

On *GPU* (NVIDIA Tesla K40c GPU, 12GB)
- [GeophysicalFlows.jl #b15419e](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/b15419e4fe093666c0b72cf1191328e631c5ed20): **3.79 ms** per time-step


Results for `nx=ny=1024` grid-points:

On *CPU* (2.9 GHz Intel Core i7 Macbook Pro with 16 GB 2133 MHz LPDDR3, using a single thread)
- [GeophysicalFlows.jl v0.3.0](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/v0.3.0): **39.43 ms** per time-step
- [pyqg v0.2.0](https://github.com/pyqg/pyqg/tree/v0.2.0): **55.76 ms** per time-step

On *GPU* (NVIDIA Tesla K40c GPU, 12GB)
- [GeophysicalFlows.jl #b15419e](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/b15419e4fe093666c0b72cf1191328e631c5ed20): **14.79 ms** per time-step
