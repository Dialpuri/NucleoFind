project(clipper-ccp4)

add_library(clipper-ccp4 STATIC 
${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_map_io.cpp
${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_mtz_io.cpp
${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_mtz_types.cpp
${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_utils.cpp
)

target_include_directories(clipper-ccp4 PRIVATE ${WRK_DIR}/checkout/libccp4 ${WRK_DIR}/checkout/clipper ${WRK_DIR}/fftw ${WRK_DIR}/rfftw)
