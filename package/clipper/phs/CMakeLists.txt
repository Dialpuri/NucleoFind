project(clipper-phs)

add_library(clipper-phs STATIC 
${WRK_DIR}/checkout/clipper/clipper/phs/phs_io.cpp
)

target_include_directories(clipper-phs PRIVATE ../.. ${WRK_DIR}/checkout/clipper ${WRK_DIR}/fftw ${WRK_DIR}/rfftw)

set(clipper-phs_HEADERS
${WRK_DIR}/checkout/clipper//clipper/phs/phs_io.h
)

target_compile_options(clipper-phs PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-phs PROPERTIES PUBLIC_HEADER "${clipper-phs_HEADERS}")

install(TARGETS clipper-phs
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/phs
)
