cmake_minimum_required(VERSION 3.0)
project(ABC_SMC_Algorithm)
set(CMAKE_CXX_STANDARD 20)

option(ABC_SMC_Algorithm_MATLABPLOTCPP "Build ${PROJECT_NAME} and the projects it depends on." OFF)
option(ABC_SMC_Algorithm_JSONCONS "Build ${PROJECT_NAME} and the projects it depends on." ON)

if(ABC_SMC_Algorithm_MATLABPLOTCPP)
    include("${CMAKE_CURRENT_SOURCE_DIR}/include/third-party-library/matplotlib-cpp/CMakeLists.txt")
endif()

#if(ABC_SMC_Algorithm_JSONCONS)
#    include("${CMAKE_CURRENT_SOURCE_DIR}/include/third-party-library/jsoncons/CMakeLists.txt")
#endif()

add_executable(ABC_SMC_Algorithm src/main.cpp 
        src/util.h src/util.cpp
        src/probability_distribution.h src/probability_distribution.cpp
        src/models/simulator.h src/models/simulator.cpp
        src/models/simulatorparameters.h src/models/simulatorparameters.cpp
        src/abc_smc_fit.cpp src/abc_smc_fit.h 
        src/models/spectrum_model.cpp src/models/spectrum_model.h src/models/spectrum.cpp
        src/pcg_random_generator.h 
        src/json_handler.cpp src/json_handler.h
        src/job_data.cpp src/config.cpp 
        src/models/example1.h src/models/example1.cpp)

# target_link_libraries(ABC_SMC_Algorithm PRIVATE matplotlib_cpp)
set_target_properties(ABC_SMC_Algorithm PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(EIGEN_DIR "./eigen-3.4.0")
find_path(EIGEN_DIR NAMES signature_of_eigen3_matrix_library
        PATHS
        /home/xxxx/Desktop/eigen-3.3.7
        PATH_SUFFIXES eigen3 eigen
        )
include_directories(${EIGEN_DIR})
