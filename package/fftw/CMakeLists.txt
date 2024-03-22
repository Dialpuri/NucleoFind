project(fftw)

add_library(fftw STATIC 
${WRK_DIR}/checkout/fftw-2.1.5/fftw/config.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/executor.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fftwf77.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fftwnd.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_1.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_10.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_11.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_12.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_13.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_14.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_15.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_16.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_2.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_3.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_32.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_4.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_5.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_6.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_64.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_7.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_8.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fn_9.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_1.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_10.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_11.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_12.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_13.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_14.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_15.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_16.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_2.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_3.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_32.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_4.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_5.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_6.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_64.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_7.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_8.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fni_9.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_10.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_16.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_2.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_3.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_32.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_4.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_5.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_6.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_64.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_7.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_8.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftw_9.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_10.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_16.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_2.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_3.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_32.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_4.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_5.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_6.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_64.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_7.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_8.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/ftwi_9.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/generic.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/malloc.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/planner.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/putils.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/rader.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/timer.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/twiddle.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/wisdom.c
${WRK_DIR}/checkout/fftw-2.1.5/fftw/wisdomio.c
)

set(fftw_HEADERS
${WRK_DIR}/checkout/fftw-2.1.5//fftw/fftw.h
${WRK_DIR}/checkout/fftw-2.1.5//fftw/fftw-int.h
${WRK_DIR}/checkout/fftw-2.1.5//fftw/config.h
)

set(fftw_HEADERS
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fftw.h
${WRK_DIR}/checkout/fftw-2.1.5/fftw/fftw-int.h
${WRK_DIR}/checkout/fftw-2.1.5/fftw/config.h
)

target_compile_options(fftw PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(fftw PROPERTIES PUBLIC_HEADER "${fftw_HEADERS}")

install(TARGETS fftw
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/fftw
)
