project(clipper-minimol)

add_library(clipper-minimol STATIC
        ${WRK_DIR}/checkout/clipper/clipper/minimol/container_minimol.cpp
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol.cpp
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_data.cpp
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_io_gemmi.cpp
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_io_mmdb.cpp
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_io_seq.cpp
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_seq.cpp
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_utils.cpp
)

target_include_directories(clipper-minimol PRIVATE
        ${WRK_DIR}/checkout/mmdb2/
        ${WRK_DIR}/checkout/clipper
        ${WRK_DIR}/checkout/fftw-2.1.5/fftw
        ${WRK_DIR}/checkout/fftw-2.1.5/rfftw
        ${WRK_DIR}/checkout/gemmi/include
)

set(clipper-minimol_HEADERS
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_data.h
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol.h
        ${WRK_DIR}/checkout/clipper/clipper/minimol/container_minimol.h
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_io_gemmi.h
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_io_mmdb.h
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_io_seq.h
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_seq.h
        ${WRK_DIR}/checkout/clipper/clipper/minimol/minimol_utils.h
        ${WRK_DIR}/checkout/clipper/clipper/clipper-minimol.h
)

target_compile_options(clipper-minimol PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-minimol PROPERTIES PUBLIC_HEADER "${clipper-minimol_HEADERS}")

install(TARGETS clipper-minimol
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper/minimol
)