Threading performance of GeophysicalFlows on Flatiron Popeye cluster:

N=1024  
1: 48.654 ms per time-step  
2: 47.953 ms per time-step  
4: 37.856 ms per time-step  
12: 27.95 ms per time-step  
24: 25.64 ms per time-step  
36: 34.79 ms per time-step  
48: 34.382 ms per time-step  

N=256  
1: 1.899 ms per time-step  
2: 2.964 ms per time-step  
4: 2.432 ms per time-step  
12: 1.829 ms per time-step  
24: 2.193 ms per time-step  
36: 3.35 ms per time-step  
48: 3.487 ms per time-step  

MPI performance of Dedalus on Flatiron Popeye cluster:

N=1024  
16: 14.386 ms per time-step  
24: 10.981 ms per time-step  
32: 8.783 ms per time-step  
48: 7.193 ms per time-step  

N=256  
16: 2.333 ms per time-step  
24: 2.253 ms per time-step  
32: 2.062 ms  per time-step  
48: 2.376 ms per time-step  
