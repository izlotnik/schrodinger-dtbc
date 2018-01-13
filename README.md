# Description
The provided Matlab codes allow to solve numerically the generalized time-dependent Schrödinger equation in unbounded domains. In case of variable coefficients becomes constant for large space variables, we can construct so-called _discrete transparent boundary conditions_ (DTBC) and use them to restrict used numercal schemes to a finite mesh. Both theoretical analysis of stability and the accomplished computations confirm that considered numerical methods coupled to the discrete TBC are effective even in the case of highly oscillating solutions and discontinuous potentials.

The following numerical methods are effectively implemented and coupled with discrete or semi-discrete TBC:
- a two-level symmetric in time (i.e. the Crank-Nicolson) scheme
- global Richardson extrapolation in time (high-order scheme in time)
- broad family of finite difference schemes with three-point parameter dependent averaging in space, incl.
  - no averaging
  - higher order Numerov-like scheme
  - scheme of the first order finite element method
- any order finite element in space numerical method (high-order scheme in space, only for 1D)
- splitting in potential (only for 2D)

For the detailed explanation, please refer to our publications.

## Main Publications
- A Zlotnik, I Zlotnik [Remarks on discrete and semi-discrete transparent boundary conditions for solving the time-dependent Schrödinger equation on the half-axis](https://doi.org/10.1515/rnam-2016-0005) _Russian Journal of Numerical Analysis and Mathematical Modelling_ 31 (1), 51-64 **2016**
- A Zlotnik, I Zlotnik [The high order method with discrete TBCs for solving the Cauchy problem for the 1D Schrödinger equation](https://doi.org/10.1515/cmam-2015-0007) _Computational Methods in Applied Mathematics_ 15 (2), 233-245 **2015**
- A Zlotnik, B Ducomet, I Zlotnik, A Romanova [Splitting in potential finite-difference schemes with discrete transparent boundary conditions for the time-dependent Schrödinger equation](https://doi.org/10.1007/978-3-319-10705-9_20) _Numerical Mathematics and Advanced Applications-ENUMATH 2013_, 203-211 **2015**
- B Ducomet, A Zlotnik, I Zlotnik [The splitting in potential Crank-Nicolson scheme with discrete transparent boundary conditions for the Schrödinger equation on a semi-infinite strip](https://doi.org/10.1051/m2an/2014004) _ESAIM: Mathematical Modelling and Numerical Analysis_ 48 (6), 1681-1699 **2014**
- I Zlotnik [Numerical methods for solving the generalized time-dependent Schrödinger equation in unbounded domains](https://www.researchgate.net/publication/273780246_Numerical_methods_for_solving_the_generalized_time-dependent_Schrodinger_equation_in_unbounded_domains) _PhD Thesis, Moscow Power Eng. Inst._ **2013**
- A Zlotnik, I Zlotnik [Finite element method with discrete transparent boundary conditions for the time-dependent 1D Schrödinger equation](http://dx.doi.org/10.3934/krm.2012.5.639) _Kinetic and Related Models_ 5 (3), 639-667 **2012**
- I Zlotnik [Family of finite-difference schemes with approximate transparent boundary conditions for the generalized nonstationary Schrödinger equation in a semi-infinite strip](https://doi.org/10.1134/S0965542511030122) _Computational Mathematics and Mathematical Physics_ 51 (3), 355-376 **2011**
- I Zlotnik [Computer simulation of the tunnel effect](https://www.researchgate.net/publication/273780247_Computer_simulation_of_the_tunnel_effect) _Moscow Power Engin. Inst. Bulletin_ 6, 10-28 **2010**
- I Zlotnik [On stability of a family of finite-difference schemes with approximate transparent boundary conditions for the time-dependent Schrödinger equation on the half-strip]() _Moscow Power Engin. Inst. Bulletin_, 127-144 **2009**
- B Ducomet, A Zlotnik, I Zlotnik [On a family of finite–difference schemes with approximate transparent boundary conditions for a generalized 1D Schrödinger equation](http://dx.doi.org/10.3934/krm.2009.2.151) _Kinetic and Related Models_ 2 (1), 151-179 **2009**

## Installation Steps
1. Download dependencies from [MATLAB Central](https://www.mathworks.com/matlabcentral/) (see below) and save them to the MATLAB path or to the project folder.
2. Run the script "Install.m" to update Matlab path and disable few warning messages.

## Dependencies
In order to run the code we will need to functions both available on [MATLAB Central](https://www.mathworks.com/matlabcentral/):
- sym2str by Martin Lawson https://www.mathworks.com/matlabcentral/fileexchange/19217-sym2str
- DataHash by Jan Simon https://www.mathworks.com/matlabcentral/fileexchange/31272-datahash

The first functon sym2str is widely used in the visualization functions to add values like $\theta=1/12$ to the labels on the plot. So if we do not need visualization you can omit that dependency. 

The second function DataHash is used to calculate MD5-hash value of the task parameter values while determining if the solution of that particular task has been already calculated and saved on the disk (if so it will be loaded from the disk and not to be calculated again). So far, if we do not need to save/load our solution to/from the hard drive we can omit the dependency of DataHash.

**Note**: We have been running our codes since Matlab R2007b. The latest version of the code is tested on Matlab R2017a. Please let us know if you run into any problem.

## Sample code
```Matlab

```

## More Information
For more information, please see the published papers or contact the author.

## Disclaimer
The code is free for academic/research purpose. Use it at your own risk and we are not responsible for any loss resulting from this code. Feel free to submit pull request for bug fixes.

## Author
[Ilya Zlotnik](https://scholar.google.ru/citations?hl=ru&user=gWphyBwAAAAJ), 2007 - 2014
