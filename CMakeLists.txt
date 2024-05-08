cmake_minimum_required(VERSION 3.15)
project(nblist LANGUAGES CUDA CXX)

# set c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(nblist SHARED)

# add source files
target_sources(nblist PRIVATE
    src/pywrap.cpp
    src/nbList.h
    src/nbList.cpp
    src/box.h
    src/box.cpp
    src/nbSchAlg.h
    src/nbSchAlg.cpp
    src/SchAlg/hashSchAlgGPU.cuh
    src/SchAlg/hashSchAlgGPU.cu
    src/SchAlg/hashSchAlgCPU.h
    src/SchAlg/hashSchAlgCPU.cpp
    src/SchAlg/octreeSchAlgCPU.h
    src/SchAlg/octreeSchAlgCPU.cpp 
    src/SchAlg/octreeSchAlgGPU.cuh
    src/SchAlg/octreeSchAlgGPU.cu 
    src/SchAlg/cellSchAlgCPU.h
    src/SchAlg/cellSchAlgCPU.cpp
    src/SchAlg/cellSchAlgGPU.cuh
    src/SchAlg/cellSchAlgGPU.cu
)

# Find pybind11
find_package(pybind11 REQUIRED)
if (NOT pybind11_FOUND)
    set(pybind11_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/external/pybind11-2.11.1/include/)
endif()
# add header files
target_include_directories(nblist PRIVATE
    ${pybind11_INCLUDE_DIRS}
)
add_subdirectory(tests)

# find openmp package and link against it
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nblist PUBLIC OpenMP::OpenMP_CXX)
endif()

# Find Python3 and link against it
find_package(Python3 REQUIRED COMPONENTS Development)
if (Python3_FOUND)
    target_link_libraries(nblist PUBLIC Python3::Python)
endif()

set(POSITION_INDEPENDENT_CODE ON)
set_target_properties(nblist PROPERTIES
    PREFIX ""
)
# set cuda compiler flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fopenmp -O3")
# set c++ compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")
