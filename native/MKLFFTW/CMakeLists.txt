# This project is used to wrap mkl deconvolution code
cmake_minimum_required(VERSION 3.9)
project(MKLFFTW)

find_package(OpenMP REQUIRED)

# find MKL include
FIND_PATH(MKL_INCLUDE_DIR $ENV {MKL_INCLUDE_DIR})

# find MKL lib
FIND_PATH(MKL_LIBRARY_DIR $ENV {MKL_LIBRARY_DIR})

# find Open MP library
FIND_PATH(OMP_LIBRARY_DIR $ENV {OMP_LIBRARY_DIR})

#include_directories(src ${MKL_INCLUDE_DIR})
link_directories(${MKL_LIBRARY_DIR} ${OMP_LIBRARY_DIR})

add_library(MKLFFTW src/MKLFFTW.cpp)

target_link_libraries(MKLFFTW mkl_core mkl_intel_lp64 mkl_intel_thread mkl_tbb_thread mkl_rt libiomp5md.lib)

install(TARGETS MKLFFTW DESTINATION lib)
