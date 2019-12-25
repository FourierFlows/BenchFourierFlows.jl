# BenchFourierFlows.jl
This project benchmarks FourierFlows performance against other similar packages

We use `nx=ny=512` grid-points and step forward an initial condition 500 times.

Results:

On *CPU* (2.9 GHz Intel Core i7 Macbook Pro with 16 GB 2133 MHz LPDDR3, using a single thread)
- [GeophysicalFlows.jl v0.3.0](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/v0.3.0): **5.136 sec**
- [pyqg v0.2.0](https://github.com/pyqg/pyqg/tree/v0.2.0): **6.638 sec**

On *GPU* (NVIDIA Tesla K40c GPU, 12GB)
- [GeophysicalFlows.jl v0.3.0](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/v0.3.0): **1.888 sec**
