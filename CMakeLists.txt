cmake_minimum_required(VERSION 3.15)
project(nblist LANGUAGES CUDA CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CUDA_NVCC_FLAGS 61)

# Add source files to a variable
set(SOURCES
    src/pywrap.cpp
    src/nbList.cpp
    src/box.cpp
    src/nbSchAlg.cpp
    src/SchAlg/hashSchAlgGPU.cu
    src/SchAlg/hashSchAlgCPU.cpp
    src/SchAlg/octreeSchAlgCPU.cpp 
    src/SchAlg/octreeSchAlgGPU.cu 
    src/SchAlg/cellSchAlgCPU.cpp
    src/SchAlg/cellSchAlgGPU.cu
)

# Find pybind11
find_package(pybind11 QUIET)

if (NOT pybind11_FOUND)
    add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11-2.11.1)
endif()

# Add the pybind11 module
pybind11_add_module(nblist SHARED)
# pybind11_add_module(nblist ${SOURCES})
target_sources(nblist PRIVATE ${SOURCES})

# Include directories
target_include_directories(nblist PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/SchAlg
)

# Find and link OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nblist PUBLIC OpenMP::OpenMP_CXX)
endif()

# Find and link Python3
find_package(Python3 REQUIRED COMPONENTS Development)
if (Python3_FOUND)
    target_link_libraries(nblist PUBLIC Python3::Python)
endif()

# Add subdirectory for tests
if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_library(dpnblist ${SOURCES})
    target_include_directories(dpnblist PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/SchAlg
    )
    target_link_libraries(dpnblist PUBLIC OpenMP::OpenMP_CXX Python3::Python)
    add_subdirectory(tests)
endif()

# Set position independent code
set(POSITION_INDEPENDENT_CODE ON)

# Set compiler flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fopenmp -O3 -arch=sm_${CUDA_NVCC_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")
