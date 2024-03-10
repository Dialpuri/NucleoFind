project(mmdb2)
cmake_minimum_required(VERSION 3.13)

add_library(mmdb2 STATIC 
${WRK_DIR}/checkout/mmdb2/mmdb2/hybrid_36.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_atom.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_bondmngr.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_chain.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_cifdefs.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_coormngr.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_cryst.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_ficif.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_io_file.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_io_stream.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_machine_.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_manager.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_mask.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_math_.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_math_align.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_math_bfgsmin.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_math_fft.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_math_graph.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_math_linalg.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_math_rand.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_mattype.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_mmcif_.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_model.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_root.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_rwbrook.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_selmngr.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_seqsuperpose.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_symop.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_tables.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_title.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_uddata.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_utils.cpp
${WRK_DIR}/checkout/mmdb2/mmdb2/mmdb_xml_.cpp
)

set(mmdb2_HEADERS
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_cryst.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_selmngr.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_xml_.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_manager.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_root.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_mmcif_.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_math_.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_title.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_atom.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_math_rand.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_utils.h
${WRK_DIR}/checkout/mmdb2//mmdb2/hybrid_36.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_tables.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_chain.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_math_fft.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_math_graph.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_machine_.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_uddata.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_io_stream.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_cifdefs.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_rwbrook.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_ficif.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_mattype.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_math_bfgsmin.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_model.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_math_linalg.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_coormngr.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_seqsuperpose.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_bondmngr.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_mask.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_math_align.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_symop.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_defs.h
${WRK_DIR}/checkout/mmdb2//mmdb2/mmdb_io_file.h
)
target_compile_options(mmdb2 PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(mmdb2 PROPERTIES PUBLIC_HEADER "${mmdb2_HEADERS}")

install(TARGETS mmdb2
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/mmdb2
)
