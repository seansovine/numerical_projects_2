# Numerical Projects 2

This repo will contain some additional numerical computation
projects for fun and for profit (by learning about some useful
numerical computation libraries).

## First project - Bessel functions

In this first small project we demonstrate the Boost odeint, GSL, and matplotplusplus
libraries working together. We will try the following:

- Solve the equation for a Bessel function using Boost's odeint.
- Plot the results using matplotplusplus.
- Use the Bessel function calculation from Gnu Scientific Library for comparison.
- Compare the two.

A complication of the ODE method is that the Bessel equation is singular at 0.
To work around this for the purpose of demonstration, we will start past 0 and use
the values produced by GSL as initial conditions for the ODE.

## Second project - power method with Eigen

This is a small example demonstrating the use of the Eigen matrix libary. It
perturbs a 3 x 3 diagonal matrix with unequal entries by a small random matrix,
then uses the power method to compute the perturbed first eigenvector and
prints out the results.

## Project infrastructure:

We have setup a basic CMake project. We will manage library dependencies as either
as submodules or using system packages. We have kept the build config very basic, but
will try out other options as the project expands over time.

## Dependencies

Libraries:

- matplotplusplus
- GSL
- Boost
- Eigen

Note that matplotplusplus and Eigen are included as submodules in this repo, 
but GSL and Boost are expected to be somewhere that CMake's `find_package`
can find them. On Linux they can be reliable installed from the package
manager.

## Resources

- [Stack Overflow: GSL with CMake](https://stackoverflow.com/questions/44821615/linking-gsl-in-cmakelists-txt)
- [GSL manual](https://www.gnu.org/software/gsl/doc/html/usage.html)
- [Wikipedia: Bessel functions](https://en.wikipedia.org/wiki/Bessel_function)
- [matplotplusplus docs](https://alandefreitas.github.io/matplotplusplus/plot-types/line-plots/function-plot/)

See the matplotplusplus readme file for using the library with CMake.

## Next ideas

- This was developed on Linux, but it would be nice to make it build cross-platform.
