
add_executable("dune_distrib_test" dune_distrib_test.cc)
target_link_dune_default_libraries("dune_distrib_test")

add_dune_ug_flags(dune_distrib_test)
add_dune_mpi_flags(dune_distrib_test)
