# Numerical Projects 2

This repo will contain some additional numerical computation
projects for fun and for profit (by learning good things).

## First project - Bessel functions

Our first goal here will be to put together some different numerical
libraries, and demonstrate them working together. We will try the following:

+ Solve the equation for a Bessel function using Boost's odeint.
+ Plot the results using matplotplusplus.
+ Use the Bessel function calculation from Gnu Scientific Library for comparison.
+ Compare the two

A complication of the ODE method is that the Bessel equation is singular at 0.
To work around this for the purpose of demonstration, we will start past 0 and use
the values produced by GSL as initial conditions for the ODE.

## Second project - power method with Eigen

This is a small example demonstrating the use of the Eigen matrix libary. It
perturbs a 3 x 3 diagonal matrix with unequal entries by a small random matrix,
then uses the power method to compute the perturbed first eigenvector and
prints out the results.

## Preliminary steps

Project infrastructure:

+ Setup a basic CMake project.
+ Add library dependencies as submodules or use system packages.
+ Setup build for basic examples, and expand as time goes on.

### Resources

+ [Stack Overflow: GSL with CMake](https://stackoverflow.com/questions/44821615/linking-gsl-in-cmakelists-txt)
+ [GSL manual](https://www.gnu.org/software/gsl/doc/html/usage.html)
+ [Wikipedia: Bessel functions](https://en.wikipedia.org/wiki/Bessel_function)
+ [matplotplusplus docs](https://alandefreitas.github.io/matplotplusplus/plot-types/line-plots/function-plot/)

See the matplotplusplus readme file for using the library with CMake.

## Dependencies

Libraries:

+ matplotplusplus
+ GSL
+ Boost
+ Eigen

## Next ideas

+ This was developed on Linux, but it would be nice to make it build cross-platform.
