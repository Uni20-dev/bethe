# This 512-node MPFR reference uses two Async workers. Reserve both CTest
# slots so the routine two-job CI does not oversubscribe it. ASan/UBSan takes
# several minutes even on a fast workstation; allow slower hosted runners
# without changing its precision, node count, tolerances or routine coverage.
# Use the discovered list: CTest's add_test does not populate CMake's if(TEST).
list(FIND bethe_xxx_four_spinon_thermodynamics_tests_TESTS
  FourSpinonThermodynamicsMpfr.PhysicalSliceWithIndependentKernelReference reference_index)
if(NOT reference_index EQUAL -1)
  set_tests_properties(FourSpinonThermodynamicsMpfr.PhysicalSliceWithIndependentKernelReference
    PROPERTIES PROCESSORS 2 TIMEOUT 1800)
endif()
