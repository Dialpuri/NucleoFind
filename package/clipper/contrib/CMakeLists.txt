project(clipper-contrib)

add_library(clipper-contrib STATIC 
${WRK_DIR}/checkout/clipper/clipper/contrib/convolution_search.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/edcalc.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/fffear.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/function_object_bases.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/mapfilter.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/originmatch.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/sfcalc.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/sfcalc_obs.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/sfscale.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/sfweight.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/skeleton.cpp
${WRK_DIR}/checkout/clipper/clipper/contrib/test_contrib.cpp
)

target_include_directories(clipper-contrib PRIVATE ../.. ${WRK_DIR}/checkout/clipper ${WRK_DIR}/fftw ${WRK_DIR}/rfftw)

set(clipper-contrib_HEADERS
${WRK_DIR}/checkout/clipper//clipper/contrib/originmatch.h
${WRK_DIR}/checkout/clipper//clipper/contrib/fffear.h
${WRK_DIR}/checkout/clipper//clipper/contrib/function_object_bases.h
${WRK_DIR}/checkout/clipper//clipper/contrib/test_contrib.h
${WRK_DIR}/checkout/clipper//clipper/contrib/sfscale.h
${WRK_DIR}/checkout/clipper//clipper/contrib/skeleton.h
${WRK_DIR}/checkout/clipper//clipper/contrib/sfcalc_obs.h
${WRK_DIR}/checkout/clipper//clipper/contrib/edcalc.h
${WRK_DIR}/checkout/clipper//clipper/contrib/sfweight.h
${WRK_DIR}/checkout/clipper//clipper/contrib/convolution_search.h
${WRK_DIR}/checkout/clipper//clipper/contrib/mapfilter.h
${WRK_DIR}/checkout/clipper//clipper/contrib/sfcalc.h
${WRK_DIR}/checkout/clipper//clipper/clipper-contrib.h
)

target_compile_options(clipper-contrib PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-contrib PROPERTIES PUBLIC_HEADER "${clipper-contrib_HEADERS}")

install(TARGETS clipper-contrib
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/contrib
)
