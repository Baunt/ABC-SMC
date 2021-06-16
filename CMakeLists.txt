cmake_minimum_required(VERSION 3.0)
project(ABC_SMC_Algorithm)

set(CMAKE_CXX_STANDARD 14)

add_executable(ABC_SMC_Algorithm main.cpp)

set(EIGEN_DIR "./eigen-3.3.9")
find_path(EIGEN_DIR NAMES signature_of_eigen3_matrix_library
        PATHS
        /home/xxxx/Desktop/eigen-3.3.7
        PATH_SUFFIXES eigen3 eigen
        )
include_directories(${EIGEN_DIR})