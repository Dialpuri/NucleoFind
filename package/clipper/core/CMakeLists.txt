project(clipper-core)

add_library(clipper-core STATIC 
${WRK_DIR}/checkout/clipper/clipper/core/atomsf.cpp
${WRK_DIR}/checkout/clipper/clipper/core/cell.cpp
${WRK_DIR}/checkout/clipper/clipper/core/clipper_instance.cpp
${WRK_DIR}/checkout/clipper/clipper/core/clipper_memory.cpp
${WRK_DIR}/checkout/clipper/clipper/core/clipper_message.cpp
${WRK_DIR}/checkout/clipper/clipper/core/clipper_stats.cpp
${WRK_DIR}/checkout/clipper/clipper/core/clipper_test.cpp
${WRK_DIR}/checkout/clipper/clipper/core/clipper_thread.cpp
${WRK_DIR}/checkout/clipper/clipper/core/clipper_types.cpp
${WRK_DIR}/checkout/clipper/clipper/core/clipper_util.cpp
${WRK_DIR}/checkout/clipper/clipper/core/container.cpp
${WRK_DIR}/checkout/clipper/clipper/core/container_hkl.cpp
${WRK_DIR}/checkout/clipper/clipper/core/container_map.cpp
${WRK_DIR}/checkout/clipper/clipper/core/container_types.cpp
${WRK_DIR}/checkout/clipper/clipper/core/coords.cpp
${WRK_DIR}/checkout/clipper/clipper/core/derivs.cpp
${WRK_DIR}/checkout/clipper/clipper/core/fftmap.cpp
${WRK_DIR}/checkout/clipper/clipper/core/fftmap_sparse.cpp
${WRK_DIR}/checkout/clipper/clipper/core/hkl_compute.cpp
${WRK_DIR}/checkout/clipper/clipper/core/hkl_data.cpp
${WRK_DIR}/checkout/clipper/clipper/core/hkl_datatypes.cpp
${WRK_DIR}/checkout/clipper/clipper/core/hkl_info.cpp
${WRK_DIR}/checkout/clipper/clipper/core/hkl_lookup.cpp
${WRK_DIR}/checkout/clipper/clipper/core/hkl_operators.cpp
${WRK_DIR}/checkout/clipper/clipper/core/map_interp.cpp
${WRK_DIR}/checkout/clipper/clipper/core/map_utils.cpp
${WRK_DIR}/checkout/clipper/clipper/core/nxmap.cpp
${WRK_DIR}/checkout/clipper/clipper/core/nxmap_operator.cpp
${WRK_DIR}/checkout/clipper/clipper/core/ramachandran.cpp
${WRK_DIR}/checkout/clipper/clipper/core/resol_basisfn.cpp
${WRK_DIR}/checkout/clipper/clipper/core/resol_fn.cpp
${WRK_DIR}/checkout/clipper/clipper/core/resol_targetfn.cpp
${WRK_DIR}/checkout/clipper/clipper/core/rotation.cpp
${WRK_DIR}/checkout/clipper/clipper/core/spacegroup.cpp
${WRK_DIR}/checkout/clipper/clipper/core/spacegroup_data.cpp
${WRK_DIR}/checkout/clipper/clipper/core/symop.cpp
${WRK_DIR}/checkout/clipper/clipper/core/test_core.cpp
${WRK_DIR}/checkout/clipper/clipper/core/test_data.cpp
${WRK_DIR}/checkout/clipper/clipper/core/xmap.cpp
)

target_include_directories(clipper-core PRIVATE ${WRK_DIR}/checkout/fftw-2.1.5/fftw ${WRK_DIR}/checkout/fftw-2.1.5/rfftw)
target_compile_options(clipper-core PUBLIC "-DFFTW_ENABLE_FLOAT")

set(clipper-core_HEADERS
${WRK_DIR}/checkout/clipper//clipper/core/container_hkl.h
${WRK_DIR}/checkout/clipper//clipper/core/fftmap.h
${WRK_DIR}/checkout/clipper//clipper/core/hkl_compute.h
${WRK_DIR}/checkout/clipper//clipper/core/test_data.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_test.h
${WRK_DIR}/checkout/clipper//clipper/core/coords.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_sysdep.h
${WRK_DIR}/checkout/clipper//clipper/core/nxmap_operator.h
${WRK_DIR}/checkout/clipper//clipper/core/derivs.h
${WRK_DIR}/checkout/clipper//clipper/core/spacegroup_data.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_precision.h
${WRK_DIR}/checkout/clipper//clipper/core/container_map.h
${WRK_DIR}/checkout/clipper//clipper/core/container_types.h
${WRK_DIR}/checkout/clipper//clipper/core/hkl_lookup.h
${WRK_DIR}/checkout/clipper//clipper/core/container.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_message.h
${WRK_DIR}/checkout/clipper//clipper/core/cell.h
${WRK_DIR}/checkout/clipper//clipper/core/resol_basisfn.h
${WRK_DIR}/checkout/clipper//clipper/core/hkl_data.h
${WRK_DIR}/checkout/clipper//clipper/core/test_core.h
${WRK_DIR}/checkout/clipper//clipper/core/rotation.h
${WRK_DIR}/checkout/clipper//clipper/core/hkl_operators.h
${WRK_DIR}/checkout/clipper//clipper/core/fftmap_sparse.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_types.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_memory.h
${WRK_DIR}/checkout/clipper//clipper/core/symop.h
${WRK_DIR}/checkout/clipper//clipper/core/resol_targetfn.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_stats.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_thread.h
${WRK_DIR}/checkout/clipper//clipper/core/spacegroup.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_instance.h
${WRK_DIR}/checkout/clipper//clipper/core/ramachandran.h
${WRK_DIR}/checkout/clipper//clipper/core/hkl_info.h
${WRK_DIR}/checkout/clipper//clipper/core/hkl_datatypes.h
${WRK_DIR}/checkout/clipper//clipper/core/xmap.h
${WRK_DIR}/checkout/clipper//clipper/core/nxmap.h
${WRK_DIR}/checkout/clipper//clipper/core/map_interp.h
${WRK_DIR}/checkout/clipper//clipper/core/map_utils.h
${WRK_DIR}/checkout/clipper//clipper/core/resol_fn.h
${WRK_DIR}/checkout/clipper//clipper/core/atomsf.h
${WRK_DIR}/checkout/clipper//clipper/core/clipper_util.h
)

target_compile_options(clipper-core PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-core PROPERTIES PUBLIC_HEADER "${clipper-core_HEADERS}")

install(TARGETS clipper-core
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/core
)

set(clipper-headers_HEADERS
${WRK_DIR}/checkout/clipper//clipper/clipper-cctbx.h
${WRK_DIR}/checkout/clipper//clipper/cctbx/clipper_cctbx.h
${WRK_DIR}/checkout/clipper//clipper/clipper-mmdb.h
${WRK_DIR}/checkout/clipper//clipper/clipper-cif.h
${WRK_DIR}/checkout/clipper//clipper/clipper-mmdbold.h
${WRK_DIR}/checkout/clipper//clipper/clipper-ccp4.h
${WRK_DIR}/checkout/clipper//clipper/minimal-clipper-hkl.h
${WRK_DIR}/checkout/clipper//clipper/clipper-phs.h
${WRK_DIR}/checkout/clipper//clipper/minimal-clipper-map.h
${WRK_DIR}/checkout/clipper//clipper/clipper-minimol.h
${WRK_DIR}/checkout/clipper//clipper/clipper.h
${WRK_DIR}/checkout/clipper//clipper/clipper-cns.h
${WRK_DIR}/checkout/clipper//clipper/clipper-contrib.h
${WRK_DIR}/checkout/clipper//clipper/cif/cif_data_io.h
${WRK_DIR}/checkout/clipper//clipper/phs/phs_io.h
)

install(FILES 
${clipper-headers_HEADERS}
DESTINATION include/clipper
)
