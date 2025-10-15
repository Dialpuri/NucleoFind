cmake_minimum_required(VERSION 3.14)
project(fftw3-wrapper LANGUAGES C CXX)

set(ENABLE_FLOAT ON CACHE BOOL "Build single precision FFTW" FORCE)
set(ENABLE_LONG_DOUBLE OFF CACHE BOOL "Build long double FFTW" FORCE)
set(ENABLE_THREADS ON CACHE BOOL "Enable FFTW threads" FORCE)
set(ENABLE_OPENMP OFF CACHE BOOL "Enable FFTW OpenMP" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared FFTW libs" FORCE)

set(FFTW3_SOURCE_DIR ${WRK_DIR}/checkout/fftw3)

add_subdirectory(${FFTW3_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/fftw_build)

if (TARGET fftw3f AND NOT TARGET FFTW3::fftw3f)
    add_library(FFTW3::fftw3f ALIAS fftw3f)
endif()

install(TARGETS fftw3f
        EXPORT FFTW3Targets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)

# Install headers
install(DIRECTORY ${FFTW3_SOURCE_DIR}/api/
        DESTINATION include/fftw3
        FILES_MATCHING PATTERN "*.h"
)

install(EXPORT FFTW3Targets
        NAMESPACE FFTW3::
        DESTINATION lib/cmake/fftw3
)
