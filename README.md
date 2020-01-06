# BenchFourierFlows.jl
This project benchmarks FourierFlows performance against other similar packages

<br/>
<br/>

Results for `nx=ny=256` grid-points:

On *CPU* (2.9 GHz Intel Core i7 Macbook Pro with 16 GB 2133 MHz LPDDR3, using one thread)

| package        | Time per time-step  |
| -------------  |:-------------:| 
|[GeophysicalFlows.jl v0.3.0](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/v0.3.0) | **1.58 ms** |
|[pyqg v0.3.0+1.gb777817](https://github.com/pyqg/pyqg/tree/b777817ecc34893407585894e3c7ab1a6160c19c)|  **3.10 ms** |
|matlab R2017a|  **9.28 ms** |
|[dedalus v2.1905+](https://bitbucket.org/dedalus-project/dedalus/src/e3f973ecb5d1861b54d12e4500e2298593f23e4a/)|  **22.48 ms** |


On *GPU* (NVIDIA Tesla K40c GPU, 12GB)

| package        | Time per time-step  |
| -------------  |:-------------:| 
|[GeophysicalFlows.jl #b15419e](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/b15419e4fe093666c0b72cf1191328e631c5ed20) | **1.02 ms** |
|matlab R2016b|  **0.975 ms** |

<br/>
<br/>

Results for `nx=ny=1024` grid-points:

On *CPU* (2.9 GHz Intel Core i7 Macbook Pro with 16 GB 2133 MHz LPDDR3, using one thread)

| package        | Time per time-step  |
| -------------  |:-------------:| 
|[GeophysicalFlows.jl v0.3.0](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/v0.3.0) | **44.29 ms** |
|[pyqg v0.3.0+1.gb777817](https://github.com/pyqg/pyqg/tree/b777817ecc34893407585894e3c7ab1a6160c19c)|  **68.52 ms** |
|matlab R2017a|  **271.92 ms** |
|[dedalus v2.1905+](https://bitbucket.org/dedalus-project/dedalus/src/e3f973ecb5d1861b54d12e4500e2298593f23e4a/)|  **357.77 ms** |


On *GPU* (NVIDIA Tesla K40c GPU, 12GB)

| package        | Time per time-step  |
| -------------  |:-------------:| 
| [GeophysicalFlows.jl #b15419e](https://github.com/FourierFlows/GeophysicalFlows.jl/tree/b15419e4fe093666c0b72cf1191328e631c5ed20) | **13.96 ms** |
|matlab R2016b|  **7.353 ms** |
