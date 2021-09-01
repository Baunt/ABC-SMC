# Approximate Bayesian Calculation - Sequential Monte Carlo Algorithm

## Installation
1. Required python in the system at least the 3 version. If you install python do not miss to adding into the environment variables. Also required to install matplotlib: ```pip install matplotlib```
2. Use x64 architecture if you build C++
## Documentation
Simplest C++ plotting library --> [matplotlib-cpp](https://github.com/lava/matplotlib-cpp)

##Tips and Tricks

[Measure execution time with high precision](https://www.geeksforgeeks.org/measure-execution-time-with-high-precision-in-c-c/)

[Computational complexity of mathematical operations](https://en.wikipedia.org/wiki/Computational_complexity_of_mathematical_operations)

|                |Complexity                          
|----------------|-------------------------------
|Addition        |*O(n)*
|Subtraction     |*O(n)*        
|Multiplication  |*O(n^2)*
|Division        |*O(n^2)*

- Combine constants where we can
- Instead of use simple math expression like from math library for example pow(x,2) use x*x. Integer quotient use like x*x*x*x. Not cute but effective. 
- Instead of dividing use the inverse of denominator like 
```vector[i] / constant --> vector[i] * constatn2 where constant2 = 1.0 / constatnt1```

## TODO

Check this library. Maybe it is faster than Eigen or vanilla C++
[NumPy to NumCpp](https://github.com/dpilger26/NumCpp)

Extended Doxygen: [Standardese](https://github.com/standardese/standardese)

Random number generator: [PCG](https://www.pcg-random.org/using-pcg-cpp.html)
