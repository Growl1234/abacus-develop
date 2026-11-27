abacus_disable_feature_definitions(__MLALGO)
abacus_disable_feature_definitions(__CUDA)
abacus_disable_feature_definitions(__ROCM)

abacus_copy_test_files(support)

addtest(
  TARGET MODULE_HAMILT_vdwTest
  LIBS parameter base device vdw
  SOURCES test_vdw.cpp)

if(ENABLE_DFTD4 AND ENABLE_MPI)
  add_test(
    NAME MODULE_HAMILT_vdwd4Test_2np
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:MODULE_HAMILT_vdwTest> --gtest_filter=vdwd4Test.*
      ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(ENABLE_LIBVDWXC)
  addtest(
    TARGET MODULE_HAMILT_vdw_nonloc
    LIBS planewave parameter base device libvdwxc::vdwxc
    SOURCES test_vdw_nonloc.cpp ../vdw_nonloc.cpp
            ../../module_xc/xc_grad_utils.cpp)
  set_tests_properties(MODULE_HAMILT_vdw_nonloc PROPERTIES ENVIRONMENT
                                                           "OMP_NUM_THREADS=1")
  if(VDWXC_ENABLE_MPI)
    target_compile_definitions(MODULE_HAMILT_vdw_nonloc
                               PRIVATE VDWXC_ENABLE_MPI)
    add_test(
      NAME MODULE_HAMILT_vdw_nonloc_3np
      COMMAND
        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:MODULE_HAMILT_vdw_nonloc> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      MODULE_HAMILT_vdw_nonloc_3np PROPERTIES PROCESSORS 3 ENVIRONMENT
                                              "OMP_NUM_THREADS=1")
  endif()
endif()
