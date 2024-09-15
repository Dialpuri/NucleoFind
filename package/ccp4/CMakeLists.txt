cmake_minimum_required(VERSION 3.15)
project(ccp4)
option(BUILD_SHARED_LIBS "Build shared library" ON)
option(INSTALL_DEV_FILES "Install headers, cmake configs" ON)
file(GLOB ccp4_src "${ccp4-dependencies_SOURCE_DIR}/ccp4/*.c")
# message(STATUS ${mmdb_src})
file(GLOB ccp4_headers "${ccp4-dependencies_SOURCE_DIR}/ccp4/*.h")
list(REMOVE_ITEM ccp4_src "${ccp4-dependencies_SOURCE_DIR}/ccp4/vmslibrary.c")
add_library(ccp4 STATIC ${ccp4_src} ${ccp4_headers})
# set(ccp4c_HEADERS ${ccp4_headers})

set_target_properties(ccp4 PROPERTIES PUBLIC_HEADER "${ccp4c_HEADERS}")
if(INSTALL_DEV_FILES)
  install(
    TARGETS ccp4
    LIBRARY DESTINATION lib
    PUBLIC_HEADER
      DESTINATION include/ccp4 # LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
      # PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ccp4
  )
endif()
