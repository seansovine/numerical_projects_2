# Project basic configuration

cmake_minimum_required(VERSION 3.21)

project(Numerical VERSION 1.2.0)

set(CMAKE_C_STANDARD 99)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# set(CMAKE_SOURCE_DIR src)
# set(PROJECT_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Third-party dependencies

add_subdirectory(matplotplusplus)

find_package(GSL REQUIRED)

find_package(Boost 1.74 REQUIRED)

# First example project

add_executable(first_example src/first_example/bessel.cpp)

target_link_libraries(first_example PUBLIC matplot GSL::gsl)

# Boost odeint example project

add_executable(boost_example src/boost_example/odeint.cpp)

target_link_libraries(boost_example Boost::boost)

add_executable(odeint_bessel src/boost_example/odeint_bessel.cpp)

target_link_libraries(odeint_bessel Boost::boost)
