# Project basic configuration

cmake_minimum_required(VERSION 3.21)

project(Numerical VERSION 1.2.0)

set(CMAKE_C_STANDARD 99)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Third-party dependencies

add_subdirectory(matplotplusplus)

find_package(GSL REQUIRED)

find_package(Boost 1.74 REQUIRED)

find_package (Eigen3 3.4 REQUIRED PATHS ./eigen)

# Bessel example project

add_executable(first_example src/bessel_example/bessel.cpp)

target_link_libraries(first_example PUBLIC matplot GSL::gsl Boost::boost)

# Eigen example project

add_executable(eigen_example src/eigen_example/example.cpp)

target_include_directories(eigen_example PRIVATE "eigen/")
