project(clipper-ccp4)

add_library(clipper-ccp4 STATIC 
${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_map_io.cpp
${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_mtz_io.cpp
${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_mtz_types.cpp
${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_utils.cpp
)

target_include_directories(clipper-ccp4 PRIVATE ${WRK_DIR}/checkout/libccp4
        ${WRK_DIR}/checkout/clipper ${WRK_DIR}/fftw ${WRK_DIR}/rfftw)

set(clipper-ccp4_HEADERS
    ${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_map_io.h
    ${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_mtz_io.h
    ${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_mtz_types.h
    ${WRK_DIR}/checkout/clipper/clipper/ccp4/ccp4_utils.h
)

set_target_properties(clipper-ccp4 PROPERTIES PUBLIC_HEADER "${clipper-ccp4_HEADERS}")

install(TARGETS clipper-ccp4
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper/ccp4
)

