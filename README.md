# Numerical Projects 2

This folder will contain some additional numerical computation
projects for fun and for profit (by learning good things).

## First goal

Our first goal here will be to put together some different numerical
tools. We will try the following:

+ Solve the equation for a Bessel function using Boost's odeint.
+ Plot the results using matplotplusplus.
+ Use the Bessel function calculation from Gnu Scientific Library for comparison.
+ Maybe approximate the derivative and try to verify derivative identity.

## Preliminary steps

Project infrastructure:

+ Setup a basic CMake project.
+ Add library dependencies as submodules.
+ Setup build for basic examples, and expand as time goes on.

### Resources

+ [Stack Overflow: GSL with CMake](https://stackoverflow.com/questions/44821615/linking-gsl-in-cmakelists-txt)
+ [GSL manual](https://www.gnu.org/software/gsl/doc/html/usage.html)
+ [Wikipedia: Bessel functions](https://en.wikipedia.org/wiki/Bessel_function)
+ [matplotplusplus docs](https://alandefreitas.github.io/matplotplusplus/plot-types/line-plots/function-plot/)

See the matplotplusplus readme file for using the library with CMake.

## Dependencies

Initially:

+ matplotplusplus
+ GSL
+ Boost
