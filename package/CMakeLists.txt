cmake_minimum_required(VERSION 3.15...3.26)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)

if(WIN32)
  set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_GENERATOR_PLATFORM x64)
endif()

project(nucleofind LANGUAGES CXX)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT SKBUILD)
  message(
    WARNING
      "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve -Ceditable.rebuild=true .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()

set(WRK_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

if(CCP4BUILD)
  find_package(
    CCP4
    REQUIRED
    mmdb2
    clipper-core
    clipper-mmdb
    clipper-ccp4
    clipper-cif
    clipper-minimol
    clipper-contrib
    ccp4c)
  find_package(gemmi-cpp)

  include_directories(${CCP4_INCLUDE_DIRS} ${CCP4_INCLUDE_DIRS}/clipper)

else()
  add_subdirectory(${WRK_DIR}/clipper/minimol)
  add_subdirectory(${WRK_DIR}/clipper/core)
  add_subdirectory(${WRK_DIR}/clipper/contrib)
  add_subdirectory(${WRK_DIR}/clipper/cns)
  add_subdirectory(${WRK_DIR}/clipper/mmdb)
  add_subdirectory(${WRK_DIR}/clipper/cif)
  add_subdirectory(${WRK_DIR}/clipper/ccp4)
  add_subdirectory(${WRK_DIR}/fftw)
  add_subdirectory(${WRK_DIR}/rfftw)
  add_subdirectory(${WRK_DIR}/ccp4)
  add_subdirectory(${WRK_DIR}/mmdb2)
  add_subdirectory(${WRK_DIR}/gemmi)
endif()
# message(CMAKE_BINARY_DIR="${CMAKE_BINARY_DIR}")

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR} ${WRK_DIR}/checkout/mmdb2
  ${WRK_DIR}/checkout/clipper ${WRK_DIR}/clipper
  ${WRK_DIR}/checkout/libccp4/ccp4 ${WRK_DIR}/checkout/gemmi/include
        ${WRK_DIR}/src/third-party
)

link_directories(${WRK_DIR}/mmdb2 ${WRK_DIR}/clipper ${WRK_DIR}/ccp4
                 ${WRK_DIR}/src ${CMAKE_BINARY_DIR})

add_library(
  nautilus_lib STATIC
  ${WRK_DIR}/src/cpp/nautilus-join.cpp
  ${WRK_DIR}/src/cpp/nautilus-rebuild-bases.cpp
  ${WRK_DIR}/src/cpp/nautilus-refine.cpp
  ${WRK_DIR}/src/cpp/nautilus-sequence.cpp
  ${WRK_DIR}/src/cpp/nautilus-ss-find.cpp
  ${WRK_DIR}/src/cpp/nautilus-target.cpp
  ${WRK_DIR}/src/cpp/nautilus-tidy.cpp
  ${WRK_DIR}/src/cpp/nautilus-tools.cpp
  ${WRK_DIR}/src/cpp/nautilus-util.cpp
  ${WRK_DIR}/src/cpp/nucleicacid_db.cpp
  ${WRK_DIR}/src/cpp/nucleofind/predicted-maps.cpp
  ${WRK_DIR}/src/cpp/nucleofind/nucleofind.cpp
  ${WRK_DIR}/src/cpp/nucleofind/backbone-tracing.cpp
  ${WRK_DIR}/src/cpp/nucleofind/fragment-library.cpp
  ${WRK_DIR}/src/cpp/nucleofind/input.cpp
  ${WRK_DIR}/src/cpp/nucleofind/refine.cpp
)

target_link_directories(
  nautilus_lib
  PUBLIC
  ${CMAKE_BINARY_DIR}/clipper/mmdb
  ${CMAKE_BINARY_DIR}/clipper/cif
  ${CMAKE_BINARY_DIR}/clipper/ccp4
  ${CMAKE_BINARY_DIR}/clipper/minimol
  ${CMAKE_BINARY_DIR}/clipper/core
  ${CMAKE_BINARY_DIR}/clipper/contrib
  ${CMAKE_BINARY_DIR}/fftw
  ${CMAKE_BINARY_DIR}/rfftw
  ${CMAKE_BINARY_DIR}/ccp4
  ${CMAKE_BINARY_DIR}/mmdb2
#  ${CMAKE_BINARY_DIR}/gemmi
)

target_link_libraries(
  nautilus_lib
  PUBLIC clipper-mmdb
         clipper-cif
         clipper-ccp4
         clipper-minimol
         clipper-core
         clipper-contrib
         rfftw
         fftw
         mmdb2
#         gemmi_cpp
)

if(CCP4BUILD)
  target_link_libraries(nautilus_lib PUBLIC ccp4c)
else()
  target_link_libraries(nautilus_lib PUBLIC ccp4)
endif()

# target_link_libraries( cnautilus PRIVATE nautilus_lib )

set(nautilus_HEADERS
    ${WRK_DIR}/src/cpp/nautilus-join.h
    ${WRK_DIR}/src/cpp/nautilus-rebuild-bases.h
    ${WRK_DIR}/src/cpp/nautilus-refine.h
    ${WRK_DIR}/src/cpp/nautilus-sequence.h
    ${WRK_DIR}/src/cpp/nautilus-ss-find.h
    ${WRK_DIR}/src/cpp/nautilus-target.h
    ${WRK_DIR}/src/cpp/nautilus-tidy.h
    ${WRK_DIR}/src/cpp/nautilus-tools.h
    ${WRK_DIR}/src/cpp/nautilus-util.h
    ${WRK_DIR}/src/cpp/nucleicacid_db.h
    ${WRK_DIR}/src/python/nautilus-include.h)

# Try to import all Python components potentially needed by nanobind
find_package(
  Python 3.8 REQUIRED
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

# We are now ready to compile the actual extension module
nanobind_add_module(
  nautilus_module STABLE_ABI NB_STATIC
  # Python bound functions here
  ${WRK_DIR}/src/python/nucleofind.cpp)

target_include_directories(
  nautilus_module
  PUBLIC ${WRK_DIR}/checkout/gemmi
         ${WRK_DIR}/mmdb2
         ${WRK_DIR}/clipper/mmdb
         ${WRK_DIR}/clipper/cif
         ${WRK_DIR}/clipper/core
         ${WRK_DIR}/clipper/contrib
         ${WRK_DIR}/clipper/ccp4
         ${WRK_DIR}/clipper/minimol
         ${WRK_DIR}/fftw
         ${WRK_DIR}/rfftw
         ${WRK_DIR}/ccp4)

target_link_directories(
  nautilus_module
  PUBLIC
  ${CMAKE_BINARY_DIR}/clipper/mmdb
  ${CMAKE_BINARY_DIR}/clipper/cif
  ${CMAKE_BINARY_DIR}/clipper/ccp4
  ${CMAKE_BINARY_DIR}/clipper/minimol
  ${CMAKE_BINARY_DIR}/clipper/contrib
  ${CMAKE_BINARY_DIR}/clipper/core
  ${CMAKE_BINARY_DIR}/fftw
  ${CMAKE_BINARY_DIR}/rfftw
  ${CMAKE_BINARY_DIR}/ccp4
  ${CMAKE_BINARY_DIR}/mmdb2
#  ${CMAKE_BINARY_DIR}/gemmi
)

target_link_libraries(
  nautilus_module
  PRIVATE nautilus_lib
          clipper-mmdb
          clipper-cif
          clipper-ccp4
          clipper-minimol
          clipper-contrib
          clipper-core
          rfftw
          fftw
          mmdb2
#          gemmi_cpp
)

if(CCP4BUILD)
  target_link_libraries(nautilus_module PUBLIC ccp4c)
else()
  target_link_libraries(nautilus_module PUBLIC ccp4)
endif()

# Install directive for scikit-build-core
install(TARGETS nautilus_module LIBRARY DESTINATION nucleofind)
