if(ENABLE_MPI)
  set(mpi_main ../../../source_base/test/mpi_test_main.cpp)
endif()
addtest(
  TARGET MODULE_GINT_divide_info
  LIBS parameter base device planewave
  SOURCES test_divide_info.cpp ../divide_info.cpp ../test/tmp_mocks.cpp
          ${mpi_main})
set_tests_properties(MODULE_GINT_divide_info PROPERTIES ENVIRONMENT
                                                        "OMP_NUM_THREADS=1")

if(ENABLE_MPI)
  foreach(ranks IN ITEMS 3 4)
    add_test(
      NAME MODULE_GINT_divide_info_mpi_${ranks}
      COMMAND
        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:MODULE_GINT_divide_info>
        ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      MODULE_GINT_divide_info_mpi_${ranks}
      PROPERTIES PROCESSORS ${ranks} TIMEOUT 60 ENVIRONMENT "OMP_NUM_THREADS=1")
  endforeach()
endif()
