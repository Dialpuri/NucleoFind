project(clipper-cif)

add_library(clipper-cif STATIC 
${WRK_DIR}/checkout/clipper/clipper/cif/cif_data_io.cpp
)

target_include_directories(clipper-cif PRIVATE ${WRK_DIR}/checkout/mmdb2/ ${WRK_DIR}/checkout/clipper ${WRK_DIR}/checkout/fftw-2.1.5/fftw ${WRK_DIR}/checkout/fftw-2.1.5/fftw  ${WRK_DIR}/checkout/gemmi/include)

set(clipper-cif_HEADERS
${WRK_DIR}/checkout/clipper/clipper/cif/cif_data_io.h
)

target_compile_options(clipper-cif PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-cif PROPERTIES PUBLIC_HEADER "${clipper-cif_HEADERS}")

install(TARGETS clipper-cif
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/cif
)
