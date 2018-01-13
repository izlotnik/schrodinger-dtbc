# Version History

## January, 2018
- Submitted to github repository
- Updated dependencies links to the MATLAB Central repository

## November, 2012
- Full code revision

## August, 2012
- Implemented fast convolution kernel calculation for big time step numbers

## July, 2012
- Implememted splitting in potential scheme for 2D Schrodinger equation
- Modified the main algorithm to reduce its space complexity

## March, 2010
- Implemented an option to add transparent boundary condition at the left boundary
- Modified solver function: vectors U_L and U_R do not change their size at each loop step

## February, 2010
- Optimized the implementation of the calculation of convolution in the function calcTimeLevelSolution.m 
  - Reduced the calculation time (x400 for m=2400) by switching from conv(X, Y) to conv2(X, Y, 'valid')
  - Reduced the total execution time by replacing all the calls of conj(X)' by transpose(X)
- Created optimized versions of Matlab's functions dot, dst, and idst (see calcDot, calcDST, and calcIDST respectively): halving execution time for m=2400
- Implemented a vectorized version of the fucntion getMatrix: execution time reduced in 5 times for m=400 function calls

## January, 2010
- Implemented a function to solve 2D Schrodinger equation

## September, 2009
- Implemented the ability to store the solution for only two time steps

## March, 2009
- Implemented the calculation of the convolution by the inetrnal Matlab's function conv
- Implemented LU factorisation as an option to solve the system of linear equation
- Reduced the calculation time in 4 times by implementing a function calcErrors to calculate errors on each time step
- Created a function plotErrorBehaviour to plot errors

## February, 2009
- Improved visualization functions
  - best location of a legend
  - save graphics in automatic mode
  - update font size in graphics objects

## August, 2008
- Implemented a normalized version of $R^m$
- Generalization to the case of arbitrary function $\rho(x)$

## July, 2008
- Simplified the calculation of the initial value R(0)
- Significantly reduced the execution time by introducing import/export procedure
- Added a function to compute different types of errors
- Introduced additional input parameters

## May, 2008
- Added codecs support to the multimedia files for screening the founded solution

## March, 2008
- Imroved the code for the case of fully arbitrary value of $\theta$
- Restructured code for better execution, introduced a function SOLVER

## February, 2008
- Fix the bug in the implementation of the Tridiagonal matrix algorithm for complex matrices, see the function TRISYS
- Fix the bug in calculation of the argument of a complex number by the standard Matlab's function angle(z) for complex $z\in\mathcol{C}$, see the function getArg
- Added code comments for main functions

## December, 2007
- Implemented the case of arbitrary value of $\theta$
- Improved performancy by implementing the recurrent relation for the sequence R(m) and by dropping out the function P_LEGENDRE
- Introduced the function to measure the execution time of each part of the algorithm

## October, 2007
- Initiated the project with few functions: P_LEGENDRE, R_SEQUENCE, TRISYS, B, V, and psi0.
- Created first visualization functions for 1D and 2D solutions of Schrodinger equation
- Implemented the solution of finite-difference scheme (only for $\theta=0$)
