# BenchFourierFlows.jl
This project benchmarks FourierFlows performance against other similar packages

On a 2.9 GHz Intel Core i7 Macbook Pro with 16 GB 2133 MHz LPDDR3 and using `nx=ny=512` and a single thread to step forward for 500 time-steps it takes

- GeophysicalFlows.jl v0.3.0: **5.136 sec**
- pyqg v0.2.0: **6.638 sec**
